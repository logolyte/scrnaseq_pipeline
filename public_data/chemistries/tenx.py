"""10x Chromium 3' v2, v3/v3.1, v4.

R1 = cell barcode (16 nt) + UMI + poly(dT), R2 = cDNA, I1/I2 = sample index. R1 is the read whose first 16 nt
are in a 10x whitelist; the whitelist with most matches gives the scrnaseq --protocol.
"""

import gzip
from collections import Counter
from functools import cache
from pathlib import Path

import pandas as pd

LAYOUT = {"fastq_1": "R1", "fastq_2": "R2"}
SAMPLESHEET_COLUMNS = ["protocol"]
FASTQ_NAME = "{source}"
ASK_SAMPLE = False
CELL_BARCODE_READ = "R1"

WHITELIST_DIR = Path(__file__).resolve().parents[1] / "whitelists"
# keys are scrnaseq --protocol values
PROTOCOLS = {
    "10XV2": {"whitelist": "10x_V2_barcode_whitelist.txt.gz", "umi_length": 10},
    "10XV3": {"whitelist": "10x_V3_barcode_whitelist.txt.gz", "umi_length": 12},
    "10XV4": {"whitelist": "10x_V4_barcode_whitelist.txt.gz", "umi_length": 12},
}
CB_LENGTH = 16
# sample index reads (I1, I2) are 8-10 bp; longer reads are R1 or R2
INDEX_MAX_LENGTH = 10
# 10x warns when fewer than 75% of barcodes are valid. 10x also accepts barcodes with one mismatch,
# here only exact matches count, so the same threshold is a bit stricter
MIN_WL = 0.75
# In 10x 5' kits R1 reads CB + UMI + TSO (template switch oligo) + cDNA, in 3' kits R1 reads CB + UMI + poly(dT).
# 5' kits use the v2 whitelist, so only the TSO tells 5' from 3' v2.
# A 26-28 bp R1 ends right after the UMI, so its TSO share is 0 and 5' data with such R1 cannot be recognized
TSO = "TTTCTTATATGGG"
MIN_TSO = 0.1
# the most frequent barcodes of the read sample; they belong to cells, not to empty droplets
TOP_BARCODES = 1000


@cache
def whitelists():
    """Loads the cell barcodes of each protocol on first use."""
    barcodes = {}
    for protocol, protocol_info in PROTOCOLS.items():
        with gzip.open(WHITELIST_DIR / protocol_info["whitelist"], "rt") as file:
            barcodes[protocol] = set(file.read().split())
    return barcodes


def profile(sequences):
    """Finds the 10x version whose whitelist has most read starts and returns it
    with the share of such reads, the share of reads with TSO after the UMI and the most frequent barcodes."""
    read_barcodes = [sequence[:CB_LENGTH] for sequence in sequences]
    whitelist_shares = {}
    for protocol, whitelist in whitelists().items():
        # True for each read whose barcode is in the whitelist; True counts as 1 in sum()
        is_in_whitelist = [barcode in whitelist for barcode in read_barcodes]
        whitelist_shares[protocol] = sum(is_in_whitelist) / len(sequences)
    protocol = max(whitelist_shares, key=whitelist_shares.get)

    barcode_length = CB_LENGTH + PROTOCOLS[protocol]["umi_length"]
    after_umi = [sequence[barcode_length : barcode_length + len(TSO)] for sequence in sequences]
    has_tso = [bases == TSO for bases in after_umi]
    top_barcodes = {barcode for barcode, _ in Counter(read_barcodes).most_common(TOP_BARCODES)}
    return {
        "protocol": protocol,
        "wl": round(whitelist_shares[protocol], 3),
        "tso": round(sum(has_tso) / len(sequences), 3),
        "cell_barcodes": top_barcodes,
    }


def roles(spot_files):
    """Assigns R1/R2 to files with the same molecules, or explains why they are not 10x 3'."""
    reads = spot_files[spot_files["length"] > INDEX_MAX_LENGTH]
    if len(reads) != 2:
        return None, f"{len(reads)} reads > {INDEX_MAX_LENGTH} bp"
    # R1 is the read whose starts are more often in the whitelist
    reads = reads.sort_values("wl", ascending=False)
    r1, r2 = reads.index
    r1_file = reads.loc[r1]

    # checked first: if R1 starts are not from the whitelist, the read taken as R1 is a random one,
    # and the next checks mean nothing
    if r1_file["wl"] < MIN_WL:
        return None, f"WL R1 = {r1_file['wl']} < {MIN_WL}"
    barcode_length = CB_LENGTH + PROTOCOLS[r1_file["protocol"]]["umi_length"]
    # STARsolo fills a too short R1 with N without an error
    if r1_file["length"] < barcode_length:
        return None, f"R1 {r1_file['length']} bp < CB+UMI {barcode_length} bp"
    if r1_file["tso"] >= MIN_TSO:
        return None, f"TSO after UMI = {r1_file['tso']} >= {MIN_TSO}: 5'"

    file_roles = pd.Series("I", index=spot_files.index)
    file_roles[r1] = "R1"
    file_roles[r2] = "R2"
    return file_roles, None
