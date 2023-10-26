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
    # These are listed for variants of which the effect is superceded by another variant's effect
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


def create_locus_tag_gene_dict(resistance_variants_csv):
    # Create dict from columns locus_tag and gene in resistance_variants_csv, with locus_tag as key and gene as value
    dict_locus_tag_gene = dict(
        zip(resistance_variants_csv["locus_tag"], resistance_variants_csv["gene"])
    )
    # Remove duplicates from dict
    dict_locus_tag_gene = {
        k: v for k, v in dict_locus_tag_gene.items() if v is not None
    }
    return dict_locus_tag_gene


def filter_for_resistance_genes(df_mutations, dict_locus_tag_gene):
    df_resistance_genes = df_mutations[
        df_mutations["locus_tag"].isin(dict_locus_tag_gene.keys())
    ]
    df_resistance_genes = df_resistance_genes.copy()
    df_resistance_genes["gene"] = df_resistance_genes["locus_tag"].map(
        dict_locus_tag_gene
    )
    return df_resistance_genes


def merge_resistance_genes_with_ref(df_resistance_genes, resistance_variants_csv):
    df_resistance_with_impact = df_resistance_genes.merge(
        resistance_variants_csv,
        how="left",
        left_on=["locus_tag", "gene", "ref_aa", "alt_aa"],
        right_on=["locus_tag", "gene", "ref_aa", "alt_aa"],
    )
    return df_resistance_with_impact


def rename_df_resistance_with_impact(df_resistance_with_impact, dict_rename):
    df_resistance_with_impact_renamed = df_resistance_with_impact.rename(
        columns=dict_rename
    )
    df_resistance_with_impact_renamed = df_resistance_with_impact_renamed[
        dict_rename.values()
    ]
    return df_resistance_with_impact_renamed


def filter_for_known_mutations(df_resistance_with_impact_renamed):
    # filter for non-empty impact
    df_known_mutations = df_resistance_with_impact_renamed[
        df_resistance_with_impact_renamed["impact"].notnull()
    ]
    return df_known_mutations


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

    locus_tag_gene_dict = create_locus_tag_gene_dict(resistance_variants_csv)

    df_resistance_genes = filter_for_resistance_genes(df_mutations, locus_tag_gene_dict)

    df_resistance_with_impact = merge_resistance_genes_with_ref(
        df_resistance_genes, resistance_variants_csv
    )

    df_all_mutations_resistance_genes = rename_df_resistance_with_impact(
        df_resistance_with_impact, dict_col_rename
    )

    df_known_resistance_mutations = filter_for_known_mutations(
        df_all_mutations_resistance_genes
    )

    df_known_resistance_mutations.to_csv(args.output, sep="\t", index=False)
    df_all_mutations_resistance_genes.to_csv(args.full_output, sep="\t", index=False)


if __name__ == "__main__":
    main()
