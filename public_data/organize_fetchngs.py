#!/usr/bin/env python3
"""Moves fetchngs FASTQs into DATASET/<GSE>/ and writes a samplesheet for nf-core/scrnaseq
with the BioSample attributes of each sample as sample_* columns.

For each GSM, shows the files with read roles from each chemistry module and asks which files to take.

A chemistry module in chemistries/ tells this script:
- LAYOUT: which read goes to which samplesheet column, e.g. {"fastq_1": "R1", "fastq_2": "R2"}
- SAMPLESHEET_COLUMNS: which of its profile values go to the samplesheet, e.g. ["protocol"]
- FASTQ_NAME: how to name moved files; "{source}" keeps the fetchngs name
- ASK_SAMPLE: True if libraries from different GSMs can be one sample (Multiome GEX and ATAC)
- profile(sequences): its values for one FASTQ, e.g. the share of barcodes in a whitelist
- roles(spot_files): which file is R1, R2, ... among files with the same reads, or why none fits
"""

import argparse
import shutil
from collections import defaultdict
from itertools import batched
from pathlib import Path

import pandas as pd
from chemistries import tenx
from fastq_info import fetch_run_metadata, list_files

CHEMISTRIES = {
    "tenx": tenx,
}
TABLE_COLUMNS = [
    "file",
    "run_accession",
    "length",
    "phred30_head",
    "phred30_tail",
    "phred30_bp",
    "read_count",
    "library",
    "original_names",
]
DROPPED_COLUMNS = ["base_count", "fastq_md5", "fastq_bytes", "fastq_ftp", "fastq_galaxy", "fastq_aspera"]


def chemistry_table(files, chemistry_name, chemistry):
    """Assigns the chemistry's read roles to each group of files with the same spots, or records why they do not fit."""
    profiles = pd.DataFrame(files[chemistry_name].tolist(), index=files.index)
    files = files.join(profiles)
    # create the columns before .loc: .loc would create them with NaN instead of None
    files["role"] = None
    files["reason"] = None
    for _, spot_files in files.groupby("spots"):
        roles, reason = chemistry.roles(spot_files)
        files.loc[spot_files.index, "role"] = roles
        files.loc[spot_files.index, "reason"] = reason
    return files


def propose(files, layout):
    """Builds the file numbers to propose for one chemistry, or explains why the files do not fit it."""
    reasons = list(files["reason"].dropna().unique())
    if reasons:
        return None, reasons
    # runs without a library ID cannot be grouped: top-up runs of one library would become separate samples
    without_library_id = files["library"] == files["run_accession"]
    if files["spots"].nunique() > 1 and without_library_id.any():
        return None, ["several runs without a library ID in original names"]

    # file numbers of each library in LAYOUT order, e.g. {"LIB1": [1, 2, 3, 4]} for R1 R2 of two runs
    file_numbers_by_library = defaultdict(list)
    for _, spot_files in files.groupby("spots", sort=False):
        library = spot_files["library"].iloc[0]
        for role in layout.values():
            file_number = spot_files.index[spot_files["role"] == role][0]
            file_numbers_by_library[library].append(file_number)

    answers_per_library = []
    for file_numbers in file_numbers_by_library.values():
        answers_per_library.append(" ".join(map(str, file_numbers)))
    return ", ".join(answers_per_library), []


def ask(files, files_by_chemistry):
    """Shows the files and each chemistry's proposal and returns the user's answer."""
    print(files[TABLE_COLUMNS].to_string())

    proposals = []
    for chemistry_name, chemistry in CHEMISTRIES.items():
        chemistry_files = files_by_chemistry[chemistry_name]
        # the names of the chemistry's profile values, e.g. ["protocol", "wl", "tso"]
        profile_names = list(files[chemistry_name].iloc[0].keys())
        roles = " ".join(chemistry.LAYOUT.values())
        print(f"\n{chemistry_name} ({roles}):")
        print(chemistry_files[["role"] + profile_names].to_string())
        file_numbers, reasons = propose(chemistry_files, chemistry.LAYOUT)
        if file_numbers is None:
            print("  not proposed: " + "; ".join(reasons))
        else:
            proposals.append(f"{chemistry_name} {file_numbers}")

    # Enter takes a proposal only if one chemistry fits; with several the user types one of them
    if len(proposals) > 1:
        print("several chemistries fit: " + "; ".join(proposals))
    proposal = proposals[0] if len(proposals) == 1 else "-"
    prompt = (
        "File numbers in role order: 'tenx 1 2' = one library; 'tenx 1 2 3 4' = one library, two runs;\n"
        "'tenx 1 2, 3 4' = two libraries; '-' = skip this GSM.\n"
        f"Enter = {proposal}\n> "
    )
    return input(prompt).strip() or proposal


def parse_answer(answer):
    """Splits the answer into chemistry and libraries of units (one unit = files of one samplesheet row):
    'tenx 1 2 3 4, 5 6' -> ('tenx', [[(1, 2), (3, 4)], [(5, 6)]])."""
    chemistry_name, numbers_text = answer.split(maxsplit=1)
    unit_size = len(CHEMISTRIES[chemistry_name].LAYOUT)
    libraries = []
    for library_text in numbers_text.split(","):
        file_numbers = [int(number) for number in library_text.split()]
        libraries.append(list(batched(file_numbers, unit_size, strict=True)))
    return chemistry_name, libraries


def library_names(gsm, libraries, files):
    """Names each library GSM, or GSM_<library> if there are several."""
    if len(libraries) == 1:
        return [gsm]
    names = []
    for units in libraries:
        first_file_number = units[0][0]
        names.append(f"{gsm}_{files.loc[first_file_number, 'library']}")
    # scrnaseq merges rows with the same sample into one run, so two libraries need different names
    if len(set(names)) != len(names):
        raise ValueError(f"{gsm}: libraries share a name {names}; group their units together")
    return names


def unit_rows(gsm, chemistry, libraries, files, samplesheet_columns, gse_folder):
    """Builds a samplesheet row for each unit and lists where its files move."""
    names = library_names(gsm, libraries, files)
    samplesheet_rows = []
    file_moves = []
    for library_name, units in zip(names, libraries, strict=True):
        sample = library_name
        if chemistry.ASK_SAMPLE:
            prompt = f"Sample for {library_name}; type the same for libraries of one sample. Enter = {library_name}\n> "
            sample = input(prompt).strip() or library_name
        for lane, unit in enumerate(units, start=1):
            # a list, because .loc reads a tuple as (row, column)
            unit_files = files.loc[list(unit)]
            # the files the user typed for one row must come from the same reads
            if unit_files["read_count"].nunique() > 1:
                raise ValueError(f"{list(unit_files['file'])} have different read counts")
            first_file = unit_files.iloc[0]
            row = first_file[samplesheet_columns + chemistry.SAMPLESHEET_COLUMNS].to_dict()
            # SRA marks R1 and R2 uploaded as separate runs as SINGLE
            row["library_layout"] = "PAIRED"
            row["sample"] = sample
            # the user typed the files in LAYOUT order, e.g. R1 then R2
            for (column, role), source in zip(chemistry.LAYOUT.items(), unit_files["path"]):
                file_name = chemistry.FASTQ_NAME.format(
                    source=source.name, prefix=library_name, lane=lane, role=role
                )
                destination = gse_folder / file_name
                row[column] = str(destination)
                file_moves.append((source, destination))
            samplesheet_rows.append(row)
    return samplesheet_rows, file_moves


def organize_gsm(gse, gsm, gsm_files, samplesheet_columns, gse_folder):
    """Asks about one GSM and returns its samplesheet rows and file moves."""
    print(f"\n{gse} {gsm}")
    files = gsm_files.reset_index(drop=True)
    # number the files from 1: the user types these numbers
    files.index += 1
    files_by_chemistry = {}
    for chemistry_name, chemistry in CHEMISTRIES.items():
        files_by_chemistry[chemistry_name] = chemistry_table(files, chemistry_name, chemistry)
    answer = ask(files, files_by_chemistry)
    if answer == "-":
        return [], []
    chemistry_name, libraries = parse_answer(answer)
    chemistry = CHEMISTRIES[chemistry_name]
    chemistry_files = files_by_chemistry[chemistry_name]
    return unit_rows(gsm, chemistry, libraries, chemistry_files, samplesheet_columns, gse_folder)


def organize_gse(gse, gse_runs, fastq_dir, samplesheet_columns, gse_folder):
    """Collects the samplesheet rows and file moves of all GSMs of one GSE."""
    original_names, attributes_by_run = fetch_run_metadata(list(gse_runs["run_accession"]))
    files = list_files(fastq_dir, gse_runs, original_names, CHEMISTRIES)
    samplesheet_rows = []
    file_moves = []
    for gsm, gsm_files in files.groupby("sample_alias"):
        gsm_rows, gsm_moves = organize_gsm(gse, gsm, gsm_files, samplesheet_columns, gse_folder)
        samplesheet_rows.extend(gsm_rows)
        file_moves.extend(gsm_moves)
    for row in samplesheet_rows:
        # BioSample attributes like genotype become columns like sample_genotype, next to sample_title of fetchngs
        for tag, value in attributes_by_run[row["run_accession"]].items():
            row[f"sample_{tag}"] = value
    return samplesheet_rows, file_moves


def move_files(file_moves, samplesheet_rows, gse_folder):
    """Moves the files and writes samplesheet.csv."""
    gse_folder.mkdir(parents=True, exist_ok=True)
    for source, destination in file_moves:
        shutil.move(source, destination)
    pd.DataFrame(samplesheet_rows).to_csv(gse_folder / "samplesheet.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetchngs", type=Path, required=True, help="fetchngs --outdir")
    parser.add_argument("--destination", type=Path, required=True, help="Destination root")
    args = parser.parse_args()

    fastq_dir = args.fetchngs / "fastq"
    fetchngs_samplesheet = pd.read_csv(
        args.fetchngs / "samplesheet" / "samplesheet.csv", dtype=str, keep_default_na=False
    )
    samplesheet_columns = list(fetchngs_samplesheet.columns.drop(DROPPED_COLUMNS))
    # we need absolute paths in the samplesheet
    destination = args.destination.resolve()

    for gse, gse_runs in fetchngs_samplesheet.groupby("study_alias"):
        gse_folder = destination / gse
        if (gse_folder / "samplesheet.csv").exists():
            print(f"{gse}: already organized, skipped")
            continue
        samplesheet_rows, file_moves = organize_gse(gse, gse_runs, fastq_dir, samplesheet_columns, gse_folder)
        if not samplesheet_rows:
            print(f"{gse}: nothing to move")
            continue
        move_files(file_moves, samplesheet_rows, gse_folder)


if __name__ == "__main__":
    main()
