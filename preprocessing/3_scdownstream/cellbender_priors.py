#!/usr/bin/env python3
"""Chooses CellBender's --expected-cells and --total-droplets-included for one sample
(priors_for_sample, used by prepare_scdownstream.py).

    expected_cells = 0.75 x the cells EmptyDrops called
    total_droplets = the barcodes with at least as many UMI as the empty-droplet peak,
                     the tallest peak of the log10-UMI histogram of the non-cell barcodes with >= 10 UMI;
                     EmptyDrops cells + 20000 if there is no usable peak

Warning codes of the 'warnings' column (details go to stderr):
  FALLBACK_NO_PEAK           no peak, or a peak below 20 UMI; used EmptyDrops cells + 20000
  FALLBACK_PEAK_ABOVE_CELLS  the peak is above the EmptyDrops cells; same fallback
  EMPTYDROPS_CALLS_EMPTIES   the lowest expected cell has < 2x the empty prior: EmptyDrops called empty droplets
  POOR_SEPARATION            cell prior < 3x empty prior: cells and empty droplets overlap; inspect the sample
  SECOND_EMPTY_PEAK          another non-cell peak is at least half as tall; the empty prior may be on the wrong one
  LARGE_TOTAL                total_droplets > EmptyDrops cells + 60000: CellBender is slow and may run out of GPU memory
"""

import sys

import anndata
import h5py
import numpy as np

EXPECTED_FRACTION = 0.75
MIN_UMI = 10  # barcodes below this are ignored when looking for the empty peak
MIN_PEAK_UMI = 20  # a peak below this is too close to MIN_UMI to trust
BIN_WIDTH = 0.05  # log10 units
SMOOTH_SIGMA = 2  # bins
FALLBACK_OFFSET = 20000
EMPTY_PRIOR_WINDOW = 500  # CellBender's empty prior is the median of this many barcodes above total_droplets
MIN_CELL_TO_EMPTY = 3
MIN_LAST_CELL_TO_EMPTY = 2
SECOND_PEAK_HEIGHT = 0.5  # relative to the tallest peak
SECOND_PEAK_DISTANCE = 0.3  # log10 units from the tallest peak
LARGE_TOTAL_OFFSET = 60000
CHUNK = 50_000_000  # matrix entries read at a time


def read_barcodes(path):
    """Returns the barcodes (obs index) of an h5ad without loading its counts."""
    with h5py.File(path, "r") as f:
        return anndata.io.read_elem(f["obs"]).index


def umi_per_barcode(path):
    """Returns the total counts of every barcode of an h5ad, reading X in chunks."""
    with h5py.File(path, "r") as f:
        X = f["X"]
        encoding = X.attrs["encoding-type"]
        # scrnaseq writes csc_matrix, where 'indices' holds the barcode of each count
        if encoding != "csc_matrix":
            raise ValueError(f"{path}: X has encoding {encoding}, expected csc_matrix")
        n_barcodes = X.attrs["shape"][0]
        data, indices = X["data"], X["indices"]
        umi = np.zeros(n_barcodes)
        for start in range(0, len(data), CHUNK):
            end = start + CHUNK
            umi += np.bincount(indices[start:end], weights=data[start:end], minlength=n_barcodes)
    return umi


def empty_peaks(umi):
    """Returns the peaks of the smoothed log10-UMI histogram as (UMI, height), tallest first."""
    x = np.log10(umi[umi >= MIN_UMI])
    if len(x) == 0:
        return []
    edges = np.arange(np.log10(MIN_UMI), x.max() + 2 * BIN_WIDTH, BIN_WIDTH)
    counts, _ = np.histogram(x, edges)
    offsets = np.arange(-3 * SMOOTH_SIGMA, 3 * SMOOTH_SIGMA + 1)
    kernel = np.exp(-0.5 * (offsets / SMOOTH_SIGMA) ** 2)
    density = np.convolve(counts, kernel / kernel.sum(), mode="same")
    peaks = [
        i for i in range(1, len(density) - 1) if density[i] > density[i - 1] and density[i] >= density[i + 1]
    ]
    centres = (edges[:-1] + edges[1:]) / 2
    return sorted(((10 ** centres[i], density[i]) for i in peaks), key=lambda p: -p[1])


def choose_priors(sample, barcodes, umi, called):
    """Chooses the settings of one sample and returns them as one row of cellbender_priors.tsv."""
    warnings = []

    def warn(code, message):
        warnings.append(code)
        print(f"{sample}: WARNING {code}: {message}", file=sys.stderr)

    is_called = barcodes.isin(called)
    n_emptydrops = len(called)
    if n_emptydrops == 0:
        raise ValueError(f"{sample}: EmptyDrops called no cells")

    expected_cells = int(round(EXPECTED_FRACTION * n_emptydrops))
    sorted_umi = np.sort(umi)[::-1]

    peaks = empty_peaks(umi[~is_called])
    peak = peaks[0][0] if peaks else None
    method = "peak"
    if peak is None or peak < MIN_PEAK_UMI:
        found = "none found" if peak is None else f"tallest at {peak:.0f} UMI"
        warn(
            "FALLBACK_NO_PEAK",
            f"no usable empty-droplet peak ({found}); total_droplets = EmptyDrops cells + {FALLBACK_OFFSET}",
        )
        method = "fallback"
    else:
        total_droplets = int((umi >= peak).sum())
        if total_droplets <= n_emptydrops:
            warn(
                "FALLBACK_PEAK_ABOVE_CELLS",
                f"the empty peak ({peak:.0f} UMI) is above the EmptyDrops cells; "
                f"total_droplets = EmptyDrops cells + {FALLBACK_OFFSET}",
            )
            method = "fallback"
        others = []
        for other_umi, height in peaks[1:]:
            is_tall = height >= SECOND_PEAK_HEIGHT * peaks[0][1]
            is_far = abs(np.log10(other_umi) - np.log10(peak)) >= SECOND_PEAK_DISTANCE
            if is_tall and is_far:
                others.append((other_umi, height))
        if others:
            listed = ", ".join(f"{u:.0f} UMI ({h / peaks[0][1]:.0%} as tall)" for u, h in others)
            warn(
                "SECOND_EMPTY_PEAK",
                f"besides the empty peak at {peak:.0f} UMI there is a peak "
                f"at {listed}; check which one is the empty droplets",
            )
    if method == "fallback":
        total_droplets = n_emptydrops + FALLBACK_OFFSET
    total_droplets = min(total_droplets, len(umi))

    if total_droplets > n_emptydrops + LARGE_TOTAL_OFFSET:
        warn(
            "LARGE_TOTAL",
            f"total_droplets {total_droplets} is more than EmptyDrops cells + "
            f"{LARGE_TOTAL_OFFSET}; CellBender will be slow and may run out of "
            f"GPU memory",
        )

    # What CellBender will compute from these (cellbender.remove_background.data.priors)
    cell_prior = float(np.exp(np.mean(np.log(np.maximum(sorted_umi[:expected_cells], 1)))))
    window = sorted_umi[max(expected_cells, total_droplets - EMPTY_PRIOR_WINDOW) : total_droplets]
    empty_prior = float(np.median(window)) if len(window) else float("nan")
    last_expected_umi = float(sorted_umi[expected_cells - 1])

    if last_expected_umi < MIN_LAST_CELL_TO_EMPTY * empty_prior:
        warn(
            "EMPTYDROPS_CALLS_EMPTIES",
            f"the lowest of the {expected_cells} expected cells has {last_expected_umi:.0f} UMI, "
            f"less than {MIN_LAST_CELL_TO_EMPTY}x the empty prior ({empty_prior:.0f}); "
            f"EmptyDrops may have called many empty droplets as cells",
        )
    if not cell_prior >= MIN_CELL_TO_EMPTY * empty_prior:
        warn(
            "POOR_SEPARATION",
            f"cell prior {cell_prior:.0f} is less than {MIN_CELL_TO_EMPTY}x "
            f"the empty prior {empty_prior:.0f}; cells and empty droplets "
            f"are poorly separated",
        )

    return {
        "sample": sample,
        "emptydrops_cells": n_emptydrops,
        "expected_cells": expected_cells,
        "total_droplets": total_droplets,
        "method": method,
        "empty_peak_umi": f"{peak:.0f}" if peak is not None else "NA",
        "cell_prior": f"{cell_prior:.0f}",
        "empty_prior": f"{empty_prior:.0f}",
        "ratio": f"{cell_prior / empty_prior:.1f}",
        "last_expected_umi": f"{last_expected_umi:.0f}",
        "warnings": ";".join(warnings) or "none",
    }


def priors_for_sample(sample, raw_path, filtered_path):
    """Returns the settings of one sample from its raw matrix and its EmptyDrops-filtered matrix."""
    barcodes = read_barcodes(raw_path)
    umi = umi_per_barcode(raw_path)
    emptydrops_barcodes = read_barcodes(filtered_path)
    return choose_priors(sample, barcodes, umi, emptydrops_barcodes)
