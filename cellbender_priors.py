#!/usr/bin/env python3
"""Choose CellBender's --expected-cells and --total-droplets-included for one sample.

    expected_cells = 0.75 x the number of cells EmptyDrops called
    total_droplets = the number of barcodes with at least as many UMI as the
                     empty-droplet peak

The empty-droplet peak is the tallest peak in a smoothed log10-UMI histogram of
the barcodes EmptyDrops did not call, using only barcodes with >= 10 UMI.
Leaving out the EmptyDrops cells keeps cells from being taken for the peak;
leaving out barcodes below 10 UMI keeps out the plateau of barcodes that only
got reads through barcode swapping. CellBender takes its empty-droplet prior
from the 500 barcodes just above total_droplets, so that prior lands on the peak.
Taking 0.75 x the EmptyDrops cells keeps borderline EmptyDrops calls out of the
cell prior; CellBender still calls more cells than expected when the data
support them.

If no usable peak is found, total_droplets falls back to EmptyDrops cells + 20000.

Warnings (written to stderr and, as short codes, to the 'warnings' column):
  FALLBACK_NO_PEAK           no peak, or a peak below 20 UMI (e.g. low-ambient
                             single-cell data); used EmptyDrops cells + 20000
  FALLBACK_PEAK_ABOVE_CELLS  the peak is above the EmptyDrops cells; same fallback
  EMPTYDROPS_CALLS_EMPTIES   the lowest expected cell has < 2x the empty prior:
                             EmptyDrops probably called many empty droplets as
                             cells, which lowers the cell prior
  POOR_SEPARATION            cell prior < 3x empty prior: cells and empty
                             droplets overlap; no prior setting fixes this, so
                             inspect the sample (damage, ambient RNA, depth)
  SECOND_EMPTY_PEAK          another peak among the non-cell barcodes is at least
                             half as tall; if the tallest is not the empty
                             droplets (e.g. debris), the empty prior is wrong
  LARGE_TOTAL                total_droplets is more than EmptyDrops cells + 60000:
                             CellBender will run much longer and may not fit in
                             GPU memory
  MISSING_BARCODES           some EmptyDrops barcodes are not in the raw matrix

Used by make_scdownstream_samplesheet.py. To check one sample by hand (needs
numpy and h5py, e.g. in the scrnaseq environment):

    python cellbender_priors.py <sample>_raw_matrix.h5ad <sample>_filtered_matrix.h5ad <sample>
"""

import sys

import h5py
import numpy as np

EXPECTED_FRACTION = 0.75
MIN_UMI = 10              # ignore barcodes below this when looking for the empty peak
MIN_PEAK_UMI = 20         # a peak below this is too close to MIN_UMI to trust
BIN_WIDTH = 0.05          # log10 units, ~12% steps in UMI
SMOOTH_SIGMA = 2          # bins
FALLBACK_OFFSET = 20000   # the previous fixed setting
EMPTY_PRIOR_WINDOW = 500  # CellBender's empty prior: median of this many barcodes above total_droplets
MIN_CELL_TO_EMPTY = 3
MIN_LAST_CELL_TO_EMPTY = 2
SECOND_PEAK_HEIGHT = 0.5  # relative to the tallest peak
SECOND_PEAK_DISTANCE = 0.3  # log10 units (2x in UMI) from the tallest peak
LARGE_TOTAL_OFFSET = 60000
CHUNK = 50_000_000        # matrix entries read at a time

COLUMNS = ["sample", "emptydrops_cells", "expected_cells", "total_droplets", "method",
           "empty_peak_umi", "cell_prior", "empty_prior", "ratio", "last_expected_umi",
           "warnings"]


def _decode(value):
    return value.decode() if isinstance(value, bytes) else value


def obs_names(f):
    """Barcodes (obs index) of an open h5ad file."""
    obs = f["obs"]
    index = obs[_decode(obs.attrs["_index"])]
    if isinstance(index, h5py.Group):  # categorical index
        return index["categories"][:][index["codes"][:]].astype(str)
    return index[:].astype(str)


def umi_per_barcode(path):
    """Barcodes and total counts per barcode of an h5ad, reading X in chunks."""
    with h5py.File(path, "r") as f:
        barcodes = obs_names(f)
        n = len(barcodes)
        X = f["X"]
        umi = np.zeros(n)
        if isinstance(X, h5py.Dataset):  # dense
            step = max(1, CHUNK // max(1, X.shape[1]))
            for start in range(0, n, step):
                umi[start:start + step] = X[start:start + step].sum(axis=1)
            return barcodes, umi
        encoding = _decode(X.attrs["encoding-type"])
        data, indices = X["data"], X["indices"]
        if encoding == "csc_matrix":  # indices are barcodes
            for start in range(0, len(data), CHUNK):
                umi += np.bincount(indices[start:start + CHUNK],
                                   weights=data[start:start + CHUNK], minlength=n)
        elif encoding == "csr_matrix":  # one row per barcode
            indptr = X["indptr"][:]
            row = 0
            while row < n:
                end = max(row + 1, int(np.searchsorted(indptr, indptr[row] + CHUNK, "right")) - 1)
                end = min(end, n)
                lo, hi = indptr[row], indptr[end]
                cumulative = np.concatenate([[0.0], np.cumsum(data[lo:hi], dtype=float)])
                umi[row:end] = cumulative[indptr[row + 1:end + 1] - lo] - cumulative[indptr[row:end] - lo]
                row = end
        else:
            raise ValueError(f"{path}: X has encoding {encoding}, expected a sparse or dense matrix")
    return barcodes, umi


def emptydrops_barcodes(path):
    with h5py.File(path, "r") as f:
        return obs_names(f)


def empty_peaks(umi):
    """Peaks of the smoothed log10-UMI histogram as (UMI, height), tallest first."""
    x = np.log10(umi[umi >= MIN_UMI])
    if len(x) == 0:
        return []
    edges = np.arange(np.log10(MIN_UMI), x.max() + 2 * BIN_WIDTH, BIN_WIDTH)
    counts, _ = np.histogram(x, edges)
    offsets = np.arange(-3 * SMOOTH_SIGMA, 3 * SMOOTH_SIGMA + 1)
    kernel = np.exp(-0.5 * (offsets / SMOOTH_SIGMA) ** 2)
    density = np.convolve(counts, kernel / kernel.sum(), mode="same")
    peaks = [i for i in range(1, len(density) - 1)
             if density[i] > density[i - 1] and density[i] >= density[i + 1]]
    centres = (edges[:-1] + edges[1:]) / 2
    return sorted(((10 ** centres[i], density[i]) for i in peaks), key=lambda p: -p[1])


def choose_priors(sample, barcodes, umi, called):
    """Settings for one sample, as a dict with the COLUMNS keys. Warnings go to stderr."""
    warnings = []

    def warn(code, message):
        warnings.append(code)
        print(f"{sample}: WARNING {code}: {message}", file=sys.stderr)

    is_called = np.isin(barcodes, called)
    n_emptydrops = len(called)
    if n_emptydrops == 0:
        raise ValueError(f"{sample}: EmptyDrops called no cells")
    if is_called.sum() != n_emptydrops:
        warn("MISSING_BARCODES", f"only {is_called.sum()} of {n_emptydrops} EmptyDrops "
                                 f"barcodes are in the raw matrix")

    expected_cells = int(round(EXPECTED_FRACTION * n_emptydrops))
    sorted_umi = np.sort(umi)[::-1]

    peaks = empty_peaks(umi[~is_called])
    peak = peaks[0][0] if peaks else None
    method = "peak"
    if peak is None or peak < MIN_PEAK_UMI:
        found = "none found" if peak is None else f"tallest at {peak:.0f} UMI"
        warn("FALLBACK_NO_PEAK", f"no usable empty-droplet peak ({found}); "
                                 f"total_droplets = EmptyDrops cells + {FALLBACK_OFFSET}")
        method = "fallback"
    else:
        total_droplets = int((umi >= peak).sum())
        if total_droplets <= n_emptydrops:
            warn("FALLBACK_PEAK_ABOVE_CELLS",
                 f"the empty peak ({peak:.0f} UMI) is above the EmptyDrops cells; "
                 f"total_droplets = EmptyDrops cells + {FALLBACK_OFFSET}")
            method = "fallback"
        others = [(u, h) for u, h in peaks[1:]
                  if h >= SECOND_PEAK_HEIGHT * peaks[0][1]
                  and abs(np.log10(u) - np.log10(peak)) >= SECOND_PEAK_DISTANCE]
        if others:
            listed = ", ".join(f"{u:.0f} UMI ({h / peaks[0][1]:.0%} as tall)" for u, h in others)
            warn("SECOND_EMPTY_PEAK", f"besides the empty peak at {peak:.0f} UMI there is a peak "
                                      f"at {listed}; check which one is the empty droplets")
    if method == "fallback":
        total_droplets = n_emptydrops + FALLBACK_OFFSET
    total_droplets = min(total_droplets, len(umi))

    if total_droplets > n_emptydrops + LARGE_TOTAL_OFFSET:
        warn("LARGE_TOTAL", f"total_droplets {total_droplets} is more than EmptyDrops cells + "
                            f"{LARGE_TOTAL_OFFSET}; CellBender will be slow and may run out of "
                            f"GPU memory")

    # What CellBender will compute from these (cellbender.remove_background.data.priors)
    cell_prior = float(np.exp(np.mean(np.log(np.maximum(sorted_umi[:expected_cells], 1)))))
    window = sorted_umi[max(expected_cells, total_droplets - EMPTY_PRIOR_WINDOW):total_droplets]
    empty_prior = float(np.median(window)) if len(window) else float("nan")
    last_expected_umi = float(sorted_umi[expected_cells - 1])

    if last_expected_umi < MIN_LAST_CELL_TO_EMPTY * empty_prior:
        warn("EMPTYDROPS_CALLS_EMPTIES",
             f"the lowest of the {expected_cells} expected cells has {last_expected_umi:.0f} UMI, "
             f"less than {MIN_LAST_CELL_TO_EMPTY}x the empty prior ({empty_prior:.0f}); "
             f"EmptyDrops may have called many empty droplets as cells")
    if not cell_prior >= MIN_CELL_TO_EMPTY * empty_prior:
        warn("POOR_SEPARATION", f"cell prior {cell_prior:.0f} is less than {MIN_CELL_TO_EMPTY}x "
                                f"the empty prior {empty_prior:.0f}; cells and empty droplets "
                                f"are poorly separated")

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
    """Settings for one sample from its raw matrix and its EmptyDrops-filtered matrix."""
    barcodes, umi = umi_per_barcode(raw_path)
    return choose_priors(sample, barcodes, umi, emptydrops_barcodes(filtered_path))


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    raw_path, filtered_path, sample = argv
    row = priors_for_sample(sample, raw_path, filtered_path)
    print("\t".join(COLUMNS))
    print("\t".join(str(row[c]) for c in COLUMNS))


if __name__ == "__main__":
    main(sys.argv[1:])
