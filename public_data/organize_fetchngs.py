#!/usr/bin/env python3
"""Moves fetchngs FASTQs into DATASET/<GSE>/ and writes a samplesheet for nf-core/scrnaseq
with the BioSample attributes of each sample as sample_* columns, and a README.txt on the source of the data.

For each GSM, shows the files with read roles from each chemistry module and asks which files to take.

A chemistry module in chemistries/ tells this script:
- LAYOUT: which read goes to which samplesheet column, e.g. {"fastq_1": "R1", "fastq_2": "R2"}
- SAMPLESHEET_COLUMNS: which of its profile values go to the samplesheet, e.g. ["protocol"]
- FASTQ_NAME: how to name moved files; "{source}" keeps the fetchngs name
- ASK_SAMPLE: True if libraries from different GSMs can be one sample (Multiome GEX and ATAC)
- CELL_BARCODE_READ: the read that starts with the cell barcode, e.g. "R1"
- profile(sequences): its values for one FASTQ, e.g. the share of barcodes in a whitelist,
  and "cell_barcodes": the most frequent cell barcodes, which tell runs of one library from different libraries
- roles(spot_files): which file is R1, R2, ... among files with the same reads, or why none fits
"""

import argparse
import datetime
import getpass
import shutil
import subprocess
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
    "original_names",
]
DROPPED_COLUMNS = ["base_count", "fastq_md5", "fastq_bytes", "fastq_ftp", "fastq_galaxy", "fastq_aspera"]
# runs sharing at least this share of their most frequent cell barcodes have the same cells.
# Runs of one library share about 80%, different libraries below 1%
SAME_CELLS = 0.5


def shared_barcodes(files, barcode_role):
    """Returns the share of cell barcodes of each run (row) found in each other run (column),
    with the file numbers of the runs' barcode reads as row and column names."""
    barcodes_of_file = files.loc[files["role"] == barcode_role, "cell_barcodes"]
    shares = pd.DataFrame(index=barcodes_of_file.index, columns=barcodes_of_file.index, dtype=float)
    for row_file, row_barcodes in barcodes_of_file.items():
        for column_file, column_barcodes in barcodes_of_file.items():
            shares.loc[row_file, column_file] = len(row_barcodes & column_barcodes) / len(row_barcodes)
    return shares.round(3)


def number_libraries(files, shares):
    """Numbers the libraries of one GSM: a run that shares most of its cell barcodes with an earlier run
    has the same cells, so it is another run of that library. Returns the library number of each file."""
    # e.g. {1: 1, 3: 1, 5: 2} for the barcode reads of two runs of library 1 and one run of library 2
    library_of_run = {}
    library_of_spots = {}
    for run_file in shares.index:
        library = 0
        for earlier_file, earlier_library in library_of_run.items():
            if shares.loc[run_file, earlier_file] >= SAME_CELLS:
                library = earlier_library
                break
        if library == 0:
            library = len(set(library_of_run.values())) + 1
        library_of_run[run_file] = library
        library_of_spots[files.loc[run_file, "spots"]] = library
    return files["spots"].map(library_of_spots)


def chemistry_table(files, chemistry_name, chemistry):
    """Assigns the chemistry's read roles to each group of files with the same spots, or records why they do not fit,
    and numbers the libraries by their cells."""
    profiles = pd.DataFrame(files[chemistry_name].tolist(), index=files.index)
    files = files.join(profiles)
    # create the columns before .loc: .loc would create them with NaN instead of None
    files["role"] = None
    files["reason"] = None
    files["library"] = None
    for _, spot_files in files.groupby("spots"):
        roles, reason = chemistry.roles(spot_files)
        files.loc[spot_files.index, "role"] = roles
        files.loc[spot_files.index, "reason"] = reason
    # libraries can be numbered only when every run has a cell barcode read
    if files["reason"].isna().all():
        shares = shared_barcodes(files, chemistry.CELL_BARCODE_READ)
        files["library"] = number_libraries(files, shares)
    return files


def propose(files, layout):
    """Builds the file numbers to propose for one chemistry, or explains why the files do not fit it."""
    reasons = list(files["reason"].dropna().unique())
    if reasons:
        return None, reasons

    # file numbers of each library in LAYOUT order, e.g. {1: [1, 2, 3, 4]} for R1 R2 of two runs
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
        # the names of the chemistry's profile values to show, e.g. ["protocol", "wl", "tso"]
        profile_names = list(files[chemistry_name].iloc[0].keys())
        profile_names.remove("cell_barcodes")
        roles = " ".join(chemistry.LAYOUT.values())
        print(f"\n{chemistry_name} ({roles}):")
        print(chemistry_files[["role"] + profile_names + ["library"]].to_string())
        file_numbers, reasons = propose(chemistry_files, chemistry.LAYOUT)
        if file_numbers is None:
            print("  not proposed: " + "; ".join(reasons))
        else:
            shares = shared_barcodes(chemistry_files, chemistry.CELL_BARCODE_READ)
            if len(shares) > 1:
                barcode_role = chemistry.CELL_BARCODE_READ
                print(f"\nShared cell barcodes of {barcode_role} files (>= {SAME_CELLS}: same library):")
                print(shares.to_string())
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


def library_names(gsm, libraries):
    """Names each library GSM, or GSM_LIB_1, GSM_LIB_2, ... in answer order if there are several."""
    if len(libraries) == 1:
        return [gsm]
    return [f"{gsm}_LIB_{number}" for number in range(1, len(libraries) + 1)]


def unit_rows(gsm, chemistry, libraries, files, samplesheet_columns, gse_folder):
    """Builds a samplesheet row for each unit and lists where its files move."""
    names = library_names(gsm, libraries)
    samplesheet_rows = []
    file_moves = []
    for library_name, units in zip(names, libraries, strict=True):
        sample = library_name
        if chemistry.ASK_SAMPLE:
            prompt = f"Sample for {library_name}; type the same for libraries of one sample. Enter = {library_name}\n> "
            sample = input(prompt).strip() or library_name
        # library numbers of the typed runs in the chemistry table
        cell_libraries = set()
        for lane, unit in enumerate(units, start=1):
            # a list, because .loc reads a tuple as (row, column)
            unit_files = files.loc[list(unit)]
            # the files the user typed for one row must come from the same reads
            if unit_files["read_count"].nunique() > 1:
                raise ValueError(f"{list(unit_files['file'])} have different read counts")
            cell_libraries.update(unit_files["library"])
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
                row[column] = file_name
                file_moves.append((source, gse_folder / file_name))
            samplesheet_rows.append(row)
        # runs of different cells in one sample would merge two channels whose cell barcodes collide
        if len(cell_libraries) > 1:
            raise ValueError(f"{library_name}: the typed runs have different cells (libraries {sorted(cell_libraries)})")
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


def write_readme(gse, gse_runs, samplesheet_rows, gse_folder):
    """Writes README.txt: the study, how its files were downloaded and sorted, by whom and when,
    and which GSMs were skipped."""
    taken_gsms = set()
    for row in samplesheet_rows:
        taken_gsms.add(row["sample_alias"])
    skipped_gsms = sorted(set(gse_runs["sample_alias"]) - taken_gsms)
    commit = subprocess.run(
        ["git", "describe", "--always", "--dirty"],
        cwd=Path(__file__).resolve().parent,
        capture_output=True,
        text=True,
        check=True,
    )
    lines = [
        f"{gse}: {gse_runs['study_title'].iloc[0]}",
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
        f"Organism: {', '.join(gse_runs['scientific_name'].unique())}",
        "",
        "Downloaded from SRA with nf-core/fetchngs (--download_method sratools)",
        f"and sorted into this folder with organize_fetchngs.py (scrnaseq_pipeline, commit {commit.stdout.strip()})",
        f"by {getpass.getuser()} on {datetime.date.today()}.",
        "",
        "samplesheet.csv   one row per sequencing run: sample (GSM), FASTQ file names, protocol,",
        "                  SRA metadata and BioSample attributes (sample_* columns)",
        "*.fastq.gz        the FASTQs chosen during sorting",
        "",
        "Index reads and skipped GSMs were left out.",
        f"Skipped GSMs: {', '.join(skipped_gsms) or 'none'}",
    ]
    (gse_folder / "README.txt").write_text("\n".join(lines) + "\n")


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

    for gse, gse_runs in fetchngs_samplesheet.groupby("study_alias"):
        gse_folder = args.destination / gse
        if (gse_folder / "samplesheet.csv").exists():
            print(f"{gse}: already organized, skipped")
            continue
        samplesheet_rows, file_moves = organize_gse(gse, gse_runs, fastq_dir, samplesheet_columns, gse_folder)
        if not samplesheet_rows:
            print(f"{gse}: nothing to move")
            continue
        move_files(file_moves, samplesheet_rows, gse_folder)
        write_readme(gse, gse_runs, samplesheet_rows, gse_folder)


if __name__ == "__main__":
    main()
