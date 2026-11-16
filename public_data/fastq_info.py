"""Describes the FASTQs of one GSE: metrics, original file names from SRA, groups of files with the same reads."""

import gzip
import os
import xml.etree.ElementTree as ET
from itertools import islice
from multiprocessing import Pool
from pathlib import Path
from statistics import median

import numpy as np
import pandas as pd
import yaml
from Bio import Entrez

# public_data/config.yaml, e.g. {"email": "you@uni.de", "skip_reads": 1000000, "sample_reads": 100000, ...}
CONFIG = yaml.safe_load((Path(__file__).resolve().parent / "config.yaml").read_text())
# sent to NCBI with each request; without it NCBI prints a warning
Entrez.email = CONFIG["email"]
SKIP_READS = CONFIG["skip_reads"]
SAMPLE_READS = CONFIG["sample_reads"]
# files read at the same time; by default one core stays free for the rest of the machine
PROCESSES = CONFIG["processes"]
if PROCESSES is None:
    PROCESSES = os.process_cpu_count() - 1
# FASTQ writes each base quality as a character whose ASCII code is Phred + 33
PHRED30_CODE = 30 + 33
# a read position counts as good when at least 70% of its bases have Phred >= 30
MIN_PHRED30_SHARE = 0.7


def fetch_run_metadata(run_accessions):
    """Fetches from SRA the original file names of each run, joined with '; ', and the BioSample attributes
    of its sample: ({"SRR1": "R1.fq.gz; R2.fq.gz"}, {"SRR1": {"genotype": "WT", "sex": "male"}})."""
    with Entrez.efetch(db="sra", id=run_accessions) as handle:
        root = ET.parse(handle).getroot()
    original_names = {}
    attributes_by_run = {}
    # one EXPERIMENT_PACKAGE holds one sample and the runs of its experiment
    for package in root.iter("EXPERIMENT_PACKAGE"):
        attributes = {}
        for attribute in package.iter("SAMPLE_ATTRIBUTE"):
            attributes[attribute.findtext("TAG")] = attribute.findtext("VALUE")
        for run_element in package.iter("RUN"):
            run = run_element.get("accession")
            file_elements = run_element.findall("SRAFiles/SRAFile[@supertype='Original']")
            file_names = [file_element.get("filename") for file_element in file_elements]
            original_names[run] = "; ".join(file_names)
            attributes_by_run[run] = attributes
    return original_names, attributes_by_run


def describe_quality(quality_lines):
    """Measures the share of Phred >= 30 bases in the first and the last quarter of the reads
    and counts the read positions where most bases have Phred >= 30."""
    # one row per read, one column per position; shorter reads are padded with code 0
    codes = np.array(quality_lines, dtype="S").view(np.uint8).reshape(len(quality_lines), -1)
    has_base = codes > 0
    is_phred30 = codes >= PHRED30_CODE
    quarter = codes.shape[1] // 4
    phred30_first_quarter = is_phred30[:, :quarter].sum() / has_base[:, :quarter].sum()
    phred30_last_quarter = is_phred30[:, -quarter:].sum() / has_base[:, -quarter:].sum()
    phred30_by_position = is_phred30.sum(axis=0) / has_base.sum(axis=0)
    return {
        "phred30_head": round(phred30_first_quarter, 2),
        "phred30_tail": round(phred30_last_quarter, 2),
        "phred30_bp": int((phred30_by_position >= MIN_PHRED30_SHARE).sum()),
    }


def describe_fastq(path):
    """Reads a sample of one FASTQ and returns its metrics and sequences."""
    first_line = 4 * SKIP_READS
    last_line = 4 * (SKIP_READS + SAMPLE_READS)
    with gzip.open(path, "rt") as file:
        lines = [line.rstrip() for line in islice(file, first_line, last_line)]
    sequences = lines[1::4]
    lengths = [len(sequence) for sequence in sequences]

    # fetchngs names files <SRX>_<SRR>_<read number>.fastq.gz, or <SRX>_<SRR>.fastq.gz for a single-file run
    run_accession = path.name.removesuffix(".fastq.gz").split("_")[1]
    row = {
        "run_accession": run_accession,
        "file": path.name,
        "path": path,
        "length": round(median(lengths)),
    }
    row.update(describe_quality(lines[3::4]))
    return row, sequences


def list_files(fastq_dir, gse_runs, original_names, chemistries):
    """Describes every FASTQ of the GSE, profiles it with each chemistry and finds files with the same reads."""
    paths = []
    for srx in gse_runs["experiment_accession"].unique():
        paths.extend(sorted(fastq_dir.glob(f"{srx}_*.fastq.gz")))
    file_rows = []
    with Pool(PROCESSES) as pool:
        for row, sequences in pool.imap(describe_fastq, paths):
            for chemistry_name, chemistry in chemistries.items():
                row[chemistry_name] = chemistry.profile(sequences)
            file_rows.append(row)
    files = pd.DataFrame(file_rows).merge(gse_runs, on="run_accession")

    # raises KeyError if NCBI returned no names for a run
    files["original_names"] = [original_names[run] for run in files["run_accession"]]

    # Some authors uploaded R1 and R2 as separate runs of one file each; such files are matched by equal read count.
    # Otherwise, R1, R2 or index reads share the same run name
    files_per_run = files.groupby("run_accession")["file"].transform("size")
    files["spots"] = np.where(files_per_run > 1, files["run_accession"], files["read_count"])
    return files
