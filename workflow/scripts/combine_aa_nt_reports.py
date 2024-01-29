#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def main(args):
    df_aa = pd.read_csv(args.aa_mutations, sep="\t")
    df_nt = pd.read_csv(args.nt_mutations, sep="\t")

    df_combined = pd.concat([df_aa, df_nt])

    df_combined.fillna("-", inplace=True)

    df_combined.sort_values(
        by=["chromosome", "position"],
        ascending=True,
        inplace=True,
    )

    df_combined.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-aa",
        "--aa-mutations",
        help="Input file with amino acid mutations",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "-nt",
        "--nt-mutations",
        help="Input file with nucleotide mutations",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file with combined mutations",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    main(args)
