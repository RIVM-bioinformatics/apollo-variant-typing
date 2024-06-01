#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Dict

import pandas as pd

dict_col_rename = {
    "gene": "genetic_element",
    "aa_change": "mutation_name",
    "impact": "impact",
    "CHROM": "chromosome",
    "POS": "position",
    "TYPE": "type_of_variant",
    "type": "type_of_consequence",
    "locus_tag": "locus_tag",
    "REF": "ref_nt",
    "ALT": "alt_nt",
    "ref_aa": "ref_aa",
    "alt_aa": "alt_aa",
    "DP": "depth",
    "AF": "allele_frequency",
}


def read_input_file(input_file: Path) -> pd.DataFrame:
    """
    Read in input file and return pandas dataframe

    Parameters
    ----------
    input_file : str
        Path to input file

    Returns
    -------
    df_input : pandas dataframe
    """
    # Open input files and filter out lines that only contain "@[0-9]+"
    # These are listed for variants of which the effect is superceded by another variant's effect
    with open(input_file, "r") as f:
        lines = f.readlines()
        lines = [
            line.rstrip("\n") for line in lines if not re.search(r"\t@[0-9]+$", line)
        ]
    # Read lines into pandas dataframe
    df_input = pd.DataFrame([line.split("\t") for line in lines[1:]])
    df_input.columns = lines[0].rstrip("\n").split("\t")
    # if AF contains a string like 0.5,0.5 convert to two rows for this record with AF 0.5
    # df_input = df_input.assign(AF=df_input["AF"].str.split(",")).explode("AF")
    # Set dtypes
    df_input = df_input.astype({"POS": int, "DP": int, "AF": float})
    df_input[["type", "locus_tag", "mutation_name"]] = df_input["BCSQ"].str.split(
        "|", expand=True
    )[[0, 1, 5]]
    df_input[["ref_aa", "alt_aa"]] = df_input["mutation_name"].str.split(
        ">", expand=True
    )
    df_input = df_input.drop(columns=["BCSQ"])
    return df_input


def create_locus_tag_gene_dict(resistance_variants_csv: pd.DataFrame) -> Dict[str, str]:
    """
    Create dictionary to map locus_tag to gene

    Parameters
    ----------
    resistance_variants_csv : pandas dataframe
        Reference CSV of AMR mutations

    Returns
    -------
    dict_locus_tag_gene : dict
        Dictionary with locus_tag as key and gene as value
    """
    # Create dict from columns locus_tag and gene in resistance_variants_csv, with locus_tag as key and gene as value
    dict_locus_tag_gene = dict(
        zip(
            resistance_variants_csv["locus_tag"],
            resistance_variants_csv["genetic_element"],
        )
    )
    # Remove duplicates from dict
    dict_locus_tag_gene = {
        k: v for k, v in dict_locus_tag_gene.items() if v is not None
    }
    return dict_locus_tag_gene


def filter_for_resistance_genes(
    df_mutations: pd.DataFrame, dict_locus_tag_gene: Dict[str, str]
) -> pd.DataFrame:
    """
    Filter df_mutations for mutations in resistance genes and add gene names

    Parameters
    ----------
    df_mutations : pandas dataframe
        Input dataframe with mutations
    dict_locus_tag_gene : dict
        Dictionary with locus_tag as key and gene as value

    Returns
    -------
    df_resistance_genes : pandas dataframe
        Dataframe with mutations in resistance genes and gene names
    """
    df_resistance_genes = df_mutations[
        df_mutations["locus_tag"].isin(dict_locus_tag_gene.keys())
    ]
    df_resistance_genes = df_resistance_genes.copy()
    df_resistance_genes["genetic_element"] = df_resistance_genes["locus_tag"].map(
        dict_locus_tag_gene
    )
    return df_resistance_genes


def merge_resistance_genes_with_ref(
    df_resistance_genes: pd.DataFrame, resistance_variants_csv: pd.DataFrame
) -> pd.DataFrame:
    """
    Add known info on resistance mutations to observed mutations

    Parameters
    ----------
    df_resistance_genes : pandas dataframe
        Dataframe with mutations in resistance genes and gene names
    resistance_variants_csv : pandas dataframe
        Reference CSV of AMR mutations

    Returns
    -------
    df_resistance_with_impact : pandas dataframe
        Dataframe with mutations in resistance genes, gene names and known info on resistance mutations
    """
    df_resistance_with_impact = df_resistance_genes.merge(
        resistance_variants_csv,
        how="left",
        left_on=["locus_tag", "genetic_element", "ref_aa", "alt_aa"],
        right_on=["locus_tag", "genetic_element", "ref_aa", "alt_aa"],
    )
    return df_resistance_with_impact


def rename_df_resistance_with_impact(
    df_resistance_with_impact: pd.DataFrame, dict_rename: Dict[str, str]
) -> pd.DataFrame:
    """
    Rename and order columns of df_resistance_with_impact

    Parameters
    ----------
    df_resistance_with_impact : pandas dataframe
        Dataframe with mutations in resistance genes, gene names and known info on resistance mutations
    dict_rename : dict
        Dictionary with old column names as keys and new column names as values

    Returns
    -------
    df_resistance_with_impact_renamed : pandas dataframe
        Dataframe with mutations in resistance genes, gene names and known info on resistance mutations, with renamed and ordered columns
    """
    df_resistance_with_impact_renamed = df_resistance_with_impact.rename(
        columns=dict_rename
    )
    df_resistance_with_impact_renamed = df_resistance_with_impact_renamed[
        dict_rename.values()
    ]
    return df_resistance_with_impact_renamed


def filter_for_known_mutations(
    df_resistance_with_impact_renamed: pd.DataFrame,
) -> pd.DataFrame:
    """
    Filter for known mutations

    Parameters
    ----------
    df_resistance_with_impact_renamed : pandas dataframe
        Dataframe with mutations in resistance genes, gene names and known info on resistance mutations, with renamed and ordered columns

    Returns
    -------
    df_known_mutations : pandas dataframe
        Dataframe with only known mutations
    """
    df_known_mutations = df_resistance_with_impact_renamed[
        df_resistance_with_impact_renamed["impact"].notnull()
    ]
    return df_known_mutations


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file", required=True, type=Path)
    parser.add_argument(
        "-o",
        "--output",
        help="Output file with only known AMR mutations",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--full-output",
        help="Output file with all mutations in resistance genes",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "-r",
        "--resistance_variants_csv",
        help="Reference CSV of AMR mutations",
        required=True,
        type=Path,
    )
    args = parser.parse_args()

    # Read in the reference list of AMR mutations
    resistance_variants_csv = pd.read_csv(args.resistance_variants_csv)

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
