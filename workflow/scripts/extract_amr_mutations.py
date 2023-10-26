#!/usr/bin/env python3

import argparse
import re

import pandas

dict_col_rename = {
    "gene": "gene",
    "aa_change": "aa_change",
    "impact": "impact",
    "CHROM": "chromosome",
    "POS": "position",
    "TYPE": "type_of_variant",
    "type": "type_of_consequence",
    "locus_tag": "locus_tag",
    "REF": "ref_nucl",
    "ALT": "alt_nucl",
    "ref_aa": "ref_aa",
    "alt_aa": "alt_aa",
    "DP": "depth",
    "AF": "allele_frequency",
}


def read_input_file(input_file):
    # Open input files and filter out lines that only contain "@[0-9]+"
    with open(input_file, "r") as f:
        lines = f.readlines()
        lines = [
            line.rstrip("\n") for line in lines if not re.match(r"^@[0-9]+$", line)
        ]
    # Read lines into pandas dataframe
    df_input = pandas.DataFrame([line.split("\t") for line in lines[1:]])
    # df_input = pandas.read_csv(input_file, sep="\t")
    df_input.columns = lines[0].rstrip("\n").split("\t")
    df_input[["type", "locus_tag", "aa_change"]] = df_input["BCSQ"].str.split(
        "|", expand=True
    )[[0, 1, 5]]
    df_input[["ref_aa", "alt_aa"]] = df_input["aa_change"].str.split(">", expand=True)
    # drop BCSQ column
    df_input = df_input.drop(columns=["BCSQ"])
    return df_input


def get_only_resistance_mutations(df_muts, df_ref, dict_rename):
    df_merged = df_muts.merge(
        df_ref,
        how="inner",
        left_on=["locus_tag", "ref_aa", "alt_aa"],
        right_on=["locus_tag", "ref_aa", "alt_aa"],
    )
    df_merged = df_merged[dict_rename.keys()]
    df_merged_renamed = df_merged.rename(columns=dict_rename)
    return df_merged_renamed


def get_allmutations_in_resistance_genes(
    df_mutations, resistance_variants_csv, dict_rename
):
    list_locus_tags = resistance_variants_csv["locus_tag"].unique()
    df_full_output = df_mutations[df_mutations["locus_tag"].isin(list_locus_tags)]
    df_full_output_rename = df_full_output.rename(columns=dict_rename)
    return df_full_output_rename


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="Output file with only known AMR mutations",
        required=True,
    )
    parser.add_argument(
        "--full-output",
        help="Output file with all mutations in resistance genes",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--resistance_variants_csv",
        help="Reference CSV of AMR mutations",
        required=True,
    )
    args = parser.parse_args()

    # Read in the reference list of AMR mutations
    resistance_variants_csv = pandas.read_csv(args.resistance_variants_csv)

    # Read in the input file
    df_mutations = read_input_file(args.input)

    df_merged = get_only_resistance_mutations(
        df_mutations, resistance_variants_csv, dict_col_rename
    )

    df_full_output = get_allmutations_in_resistance_genes(
        df_mutations, resistance_variants_csv, dict_col_rename
    )

    df_merged.to_csv(args.output, sep="\t", index=False)
    df_full_output.to_csv(args.full_output, sep="\t", index=False)


if __name__ == "__main__":
    main()
