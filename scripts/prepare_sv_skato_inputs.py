#!/usr/bin/env python3
"""Convert a chromosome TSV into typed inputs for the Terra WDL."""

import argparse
import csv
import json
from pathlib import Path


COLUMNS = ("chromosome", "pgen", "pvar", "psam", "molecular_traits")
CHROMOSOME_ALIASES = {"23": "X", "24": "Y", "26": "MT", "M": "MT"}


def read_manifest(path):
    entries = []
    seen = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if sorted(rows.fieldnames or ()) != sorted(COLUMNS):
            raise ValueError(f"Expected TSV columns: {', '.join(COLUMNS)}")
        for row in rows:
            if None in row or any(not (row.get(column) or "").strip() for column in COLUMNS):
                raise ValueError(f"Invalid manifest row at line {rows.line_num}.")
            entry = {column: row[column].strip() for column in COLUMNS}
            if any("\n" in value or "\r" in value for value in entry.values()):
                raise ValueError(f"Line breaks are invalid in manifest fields at line {rows.line_num}.")
            chromosome = entry["chromosome"].upper().removeprefix("CHR")
            chromosome = CHROMOSOME_ALIASES.get(chromosome, chromosome)
            if not chromosome or chromosome in seen:
                raise ValueError(f"Blank or repeated chromosome at line {rows.line_num}.")
            seen.add(chromosome)
            entries.append(entry)
    if not entries:
        raise ValueError("The manifest must contain at least one chromosome.")
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--cds-gtf", required=True)
    parser.add_argument("--lof-table", required=True)
    parser.add_argument("--docker-image", required=True)
    parser.add_argument("--covariates")
    parser.add_argument("--covariate-columns", help="Comma-delimited column names")
    parser.add_argument("--fdr-method", choices=("storey", "BH"), default="storey")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if bool(args.covariates) != bool(args.covariate_columns):
        parser.error("Supply --covariates and --covariate-columns together.")
    try:
        manifest = read_manifest(args.manifest)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    inputs = {
        "sv_skato.chromosome_manifest": manifest,
        "sv_skato.cds_gtf": args.cds_gtf,
        "sv_skato.lof_table": args.lof_table,
        "sv_skato.docker_image": args.docker_image,
        "sv_skato.fdr_method": args.fdr_method,
    }
    if args.covariates:
        inputs["sv_skato.covariates"] = args.covariates
        inputs["sv_skato.covariate_columns"] = args.covariate_columns
    args.output.write_text(json.dumps(inputs, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
