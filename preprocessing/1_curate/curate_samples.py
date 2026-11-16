#!/usr/bin/env python3
"""Asks which samples to analyze and which mode of config.yaml to use, e.g. sc (cells) or sn (nuclei).
With --mode it takes all samples and asks nothing.
Writes their runs to <outdir>/samplesheet.csv, a README.txt and a copy of config.yaml with the answers to <outdir>/config.yaml."""

import argparse
import datetime
import getpass
import re
import subprocess
from pathlib import Path

import pandas as pd
import yaml

CONFIG = Path(__file__).resolve().parents[1] / "config.yaml"
FASTQ_COLUMNS = ["sample", "fastq_1", "fastq_2"]


def choose_samples(samples):
    """Asks for sample numbers like '1,3,5' or '1-5' (Enter = all) and returns the chosen rows
    once the user confirms their names."""
    while True:
        answer = input("Samples (1,3,5 or 1-5; Enter = all): ").replace(" ", "")
        if answer == "":
            chosen_samples = samples
        elif re.fullmatch(r"\d+(-\d+)?(,\d+(-\d+)?)*", answer):
            numbers = set()
            for part in answer.split(","):
                start, _, end = part.partition("-")
                numbers.update(range(int(start), int(end or start) + 1))
            if not numbers <= set(samples.index):
                print("Invalid selection")
                continue
            chosen_samples = samples.loc[sorted(numbers)]
        else:
            print("Invalid selection")
            continue

        names = ", ".join(chosen_samples["sample"])
        if input(f"Take {names}? (Enter = yes, n = choose again): ").strip() == "":
            return chosen_samples


def write_readme(samplesheet_path, samples, chosen_samples, mode, protocol, outdir):
    """Writes README.txt: the study from the README.txt of public_data next to the samplesheet, if there is one,
    the chosen samples and settings, by whom and when, and what the folders hold."""
    lines = []
    study_readme = samplesheet_path.parent / "README.txt"
    if study_readme.exists():
        # the first paragraph: study title, GEO link, organism
        lines += [study_readme.read_text().split("\n\n")[0], ""]
    commit = subprocess.run(
        ["git", "describe", "--always", "--dirty"],
        cwd=Path(__file__).resolve().parent,
        capture_output=True,
        text=True,
        check=True,
    )
    skipped_samples = samples.loc[~samples["sample"].isin(chosen_samples["sample"]), "sample"]
    lines += [
        f"Preprocessed with https://github.com/logolyte/scrnaseq_pipeline (commit {commit.stdout.strip()})",
        f"by {getpass.getuser()} on {datetime.date.today()}: mode {mode}, protocol {protocol}.",
        "",
        f"Samples: {', '.join(chosen_samples['sample'])}",
        f"Skipped samples: {', '.join(skipped_samples) or 'none'}",
        "",
        "config.yaml, samplesheet.csv   settings and runs of this preprocessing",
        "slurm_logs/                    one log per step",
        "scrnaseq/, scdownstream/       Nextflow launch folders: params.yaml, .nextflow.log, the logs of each task in work/;",
        "                               scdownstream/cellbender_priors.tsv: CellBender priors and warnings",
        "results/*/pipeline_info/       trace, timeline and report of each Nextflow run",
        "results/starsolo/              FastQC, MultiQC, STAR logs and Solo.out statistics per sample",
        "results/scdownstream/          QC report, MultiQC, CellBender reports and cell barcodes, QC plots and tables",
        "results/cellassign/            cell types of all cells, marker genes",
        "velocyto/, results/velocyto*/  spliced and unspliced counts, only with run_velocyto",
        "results/final/*.h5ad           the result: log1p in X, raw and decontX counts, cell types, sample metadata",
        "",
        "The copy from the server to storage (step 5 of preprocessing/README.md) leaves out BAM, count matrices",
        "and intermediate h5/h5ad/pkl files.",
    ]
    (outdir / "README.txt").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samplesheet", type=Path, help="One row per run: sample, fastq_1, fastq_2, protocol, metadata")
    parser.add_argument("outdir", type=Path, help="Output folder of this preprocessing run")
    parser.add_argument("--mode", help="A mode of config.yaml, e.g. sn; takes all samples without asking")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    samplesheet = pd.read_csv(args.samplesheet, dtype=str, keep_default_na=False)
    for column in ["fastq_1", "fastq_2"]:
        samplesheet[column] = [str(args.samplesheet.parent.resolve() / path) for path in samplesheet[column]]
    lowercase_names = {column.lower(): column for column in samplesheet.columns}
    metadata_columns = []
    for column in config["metadata_columns"]:
        if column.lower() in lowercase_names:
            samplesheet = samplesheet.rename(columns={lowercase_names[column.lower()]: column})
            metadata_columns.append(column)
    # every sample gets one value per column: the runs of a sample must agree
    values_per_sample = samplesheet.groupby("sample")[["protocol"] + metadata_columns].nunique()
    differing = values_per_sample[(values_per_sample > 1).any(axis=1)]
    if len(differing) > 0:
        raise ValueError(f"runs of one sample have different values:\n{differing}")

    samples = samplesheet.drop_duplicates("sample")[["sample", "protocol"] + metadata_columns].reset_index(drop=True)
    samples.index += 1
    print(samples.to_string())

    if args.mode:
        chosen_samples = samples
        mode = args.mode
    else:
        chosen_samples = choose_samples(samples)
        mode = input(f"Mode ({', '.join(config['modes'])}): ").strip()
    # scrnaseq takes one --protocol per run
    if chosen_samples["protocol"].nunique() > 1:
        raise ValueError("the chosen samples have different protocols")
    if mode not in config["modes"]:
        raise ValueError(f"'{mode}' is not one of {list(config['modes'])}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    runs = samplesheet[samplesheet["sample"].isin(chosen_samples["sample"])]
    runs[FASTQ_COLUMNS + metadata_columns].to_csv(args.outdir / "samplesheet.csv", index=False)
    protocol = chosen_samples["protocol"].iloc[0]
    write_readme(args.samplesheet, samples, chosen_samples, mode, protocol, args.outdir)
    # written last: run_pipeline.py treats an existing config.yaml as a finished curation
    run_config = config | {"mode": mode, "protocol": protocol}
    (args.outdir / "config.yaml").write_text(yaml.safe_dump(run_config, sort_keys=False))


if __name__ == "__main__":
    main()
