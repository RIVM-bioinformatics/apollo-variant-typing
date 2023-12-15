#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Dict

import pandas as pd

dict_col_rename = {
    "genetic_element": "genetic_element",
    "mutation_name": "mutation_name",
    "comparison_type": "comparison_type",
    "impact": "impact",
    "CHROM": "chromosome",
    "POS": "position",
    "TYPE": "type_of_variant",
    "REF": "ref_nt",
    "ALT": "alt_nt",
    "DP": "depth",
    "AF": "allele_frequency",
}


def find_exact_matches(
    df_resistance_variants: pd.DataFrame,
    df_mutations: pd.DataFrame,
    dict_col_rename: Dict[str, str],
) -> pd.DataFrame:
    """
    Find exact matches between df_resistance_variants and df_mutations

    Parameters
    ----------
    df_resistance_variants : pandas dataframe
        Reference CSV of AMR mutations
    df_mutations : pandas dataframe
        Input dataframe with mutations

    Returns
    -------
    df_exact_matches : pandas dataframe
        Dataframe with exact matches between df_resistance_variants and df_mutations
    """
    df_exact_matches = df_mutations.merge(
        df_resistance_variants,
        how="inner",
        left_on=["CHROM", "POS", "REF", "ALT"],
        right_on=["chrom", "position", "ref_nt", "alt_nt"],
    )
    df_exact_matches_relevant_cols = df_exact_matches[dict_col_rename.keys()]
    df_exact_matches_relevant_cols_renamed = df_exact_matches_relevant_cols.rename(
        columns=dict_col_rename
    )
    return df_exact_matches_relevant_cols_renamed


def find_large_indels(
    df_mutations: pd.DataFrame, screen_region: tuple[int, int], genetic_element: str
) -> pd.DataFrame:
    """
    Find large indels in screen region

    Parameters
    ----------
    df_mutations : pandas dataframe
        Input dataframe with mutations
    screen_region : tuple[int, int]
        Screen region

    Returns
    -------
    df_large_indels : pandas dataframe
        Dataframe with large indels in screen region
    """
    df_large_indels = df_mutations[
        (df_mutations["TYPE"] == "INDEL")
        & (df_mutations["POS"] >= screen_region[0])
        & (df_mutations["POS"] <= screen_region[1])
        & (
            (df_mutations["REF"].apply(lambda x: len(x)) >= 5)
            | (df_mutations["ALT"].apply(lambda x: len(x)) >= 5)
        )
    ].copy()
    df_large_indels["genetic_element"] = genetic_element
    df_large_indels["mutation_name"] = "-"

    for index, row in df_large_indels.iterrows():
        if len(row["REF"]) > len(row["ALT"]):
            mutation_name = f"possible_tandem_repeat_length_{len(row['REF']) - 1}"
        else:
            mutation_name = f"possible_tandem_repeat_length_{len(row['ALT']) - 1}"
        df_large_indels.loc[index, "mutation_name"] = mutation_name
    return df_large_indels


def screen_for_possible_cnv_in_known_regions(
    df_resistance_variants: pd.DataFrame,
    df_mutations: pd.DataFrame,
    dict_col_rename: Dict[str, str],
) -> pd.DataFrame:
    """
    Screen for possible CNV in known regions
    """
    df_resistance_variants_tandem_repeat = df_resistance_variants[
        df_resistance_variants["comparison_type"] == "tandem_repeat"
    ]
    list_df_possible_cnvs = []
    for index, row in df_resistance_variants_tandem_repeat.iterrows():
        screen_region = (row["position"] - 50, row["position"] + 50)
        list_df_possible_cnvs.append(
            find_large_indels(df_mutations, screen_region, row["genetic_element"])
        )
    df_possible_cnvs = pd.concat(list_df_possible_cnvs).drop_duplicates()
    df_possible_cnvs["comparison_type"] = "tandem_repeat"
    df_possible_cnvs["impact"] = "unknown"

    df_possible_cnvs_relevant_cols = df_possible_cnvs[dict_col_rename.keys()]
    df_possible_cnvs_relevant_cols_renamed = df_possible_cnvs_relevant_cols.rename(
        columns=dict_col_rename
    )

    return df_possible_cnvs_relevant_cols_renamed


def combine_exact_matches_and_possible_cnvs(
    df_exact_matches: pd.DataFrame, df_possible_cnvs: pd.DataFrame
) -> pd.DataFrame:
    """
    Combine exact matches and possible CNVs

    Parameters
    ----------
    df_exact_matches : pandas dataframe
        Dataframe with exact matches between df_resistance_variants and df_mutations
    df_possible_cnvs : pandas dataframe
        Dataframe with possible CNVs in known regions

    Returns
    -------
    df_output : pandas dataframe
        Dataframe with exact matches and possible CNVs
    """
    df_possible_cnvs["sort_col"] = 0
    df_exact_matches["sort_col"] = 1

    df_combined = pd.concat([df_exact_matches, df_possible_cnvs]).sort_values(
        by=["sort_col"],
        ascending=False,
    )

    df_duplicates_removed = df_combined.drop_duplicates(
        subset=["chromosome", "position", "ref_nt", "alt_nt"],
        keep="first",
    )

    df_output = df_duplicates_removed.drop(columns=["sort_col"])
    df_output = df_output.sort_values(by=["chromosome", "position", "ref_nt", "alt_nt"])

    return df_output


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
    # In rare cases, BCSQ can be a column of only NA which will otherwise be read in as a float
    df_mutations = pd.read_csv(args.input, sep="\t", dtype={"BCSQ": object})

    df_exact_matches = find_exact_matches(
        resistance_variants_csv, df_mutations, dict_col_rename
    )

    df_possible_cnvs = screen_for_possible_cnv_in_known_regions(
        resistance_variants_csv, df_mutations, dict_col_rename
    )

    df_output = combine_exact_matches_and_possible_cnvs(
        df_exact_matches, df_possible_cnvs
    )

    df_output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
