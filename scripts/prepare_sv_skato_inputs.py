#!/usr/bin/env python3
"""Write Cromwell inputs that point to a chromosome TSV in Cloud Storage."""

import argparse
import json
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="gs:// URI of the chromosome manifest TSV")
    parser.add_argument("--cds-gtf", help="Use annotated CDS intervals for SV overlap")
    parser.add_argument("--window-padding-bp", type=int, help="Padding per side for phenotype BED intervals")
    parser.add_argument("--lof-table", required=True)
    parser.add_argument("--docker-image", help="Override the WDL's default RareSV image")
    parser.add_argument("--covariates")
    parser.add_argument("--covariate-columns", help="Comma-delimited column names")
    parser.add_argument("--fdr-method", choices=("storey", "BH"), default="storey")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if bool(args.covariates) != bool(args.covariate_columns):
        parser.error("Supply --covariates and --covariate-columns together.")
    if args.window_padding_bp is not None and not 0 <= args.window_padding_bp <= 2147483647:
        parser.error("Window padding must be an integer from 0 to 2147483647 bp.")
    if not re.fullmatch(r"gs://[^/\s]+/[^\t\r\n]+", args.manifest):
        parser.error("--manifest must be a gs:// URI for the chromosome manifest TSV.")
    inputs = {
        "sv_skato.chromosome_manifest": args.manifest,
        "sv_skato.lof_table": args.lof_table,
        "sv_skato.fdr_method": args.fdr_method,
    }
    if args.docker_image:
        inputs["sv_skato.docker_image"] = args.docker_image
    if args.cds_gtf:
        inputs["sv_skato.cds_gtf"] = args.cds_gtf
    if args.window_padding_bp is not None:
        inputs["sv_skato.window_padding_bp"] = args.window_padding_bp
    if args.covariates:
        inputs["sv_skato.covariates"] = args.covariates
        inputs["sv_skato.covariate_columns"] = args.covariate_columns
    args.output.write_text(json.dumps(inputs, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
