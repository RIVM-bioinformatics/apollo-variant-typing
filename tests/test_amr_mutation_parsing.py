import unittest
from pathlib import Path
from sys import path

import pandas as pd

from workflow.scripts.compare_aa_mutations import (
    create_locus_tag_gene_dict,
    filter_for_known_mutations,
    filter_for_resistance_genes,
    merge_resistance_genes_with_ref,
    read_input_file,
    rename_df_resistance_with_impact,
)
from workflow.scripts.compare_nt_mutations import (
    combine_exact_matches_and_possible_cnvs,
    find_exact_matches,
    find_large_indels,
    screen_for_possible_cnv_in_known_regions,
)

# init df_mutations with missense SNP and synonymous SNP in the same genetic element, and a large intergenic INDEL (promoter mutation)
df_mutations = pd.read_csv("tests/test_files/df_mutations.tsv", sep="\t")

df_aa_resistance_variants = pd.read_csv(
    "tests/test_files/df_aa_resistance_variants.tsv", sep="\t"
)

df_nt_resistance_variants = pd.read_csv(
    "tests/test_files/df_nt_resistance_variants.tsv", sep="\t"
)

dict_col_rename_nt = {
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

dict_col_rename_aa = {
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


class TestNtComparison(unittest.TestCase):
    df_exact_matches_correct = pd.read_csv(
        "tests/test_files/df_exact_matches_correct.tsv", sep="\t"
    )
    df_find_large_indels_correct = pd.read_csv(
        "tests/test_files/df_find_large_indels_correct.tsv",
        sep="\t",
        dtype={"BCSQ": object},
    )

    df_find_large_indels_no_hits_correct = pd.read_csv(
        "tests/test_files/df_find_large_indels_no_hits_correct.tsv",
        sep="\t",
        dtype={"BCSQ": object},
    )

    df_screen_for_possible_cnv_in_known_regions_correct = pd.read_csv(
        "tests/test_files/df_screen_for_possible_cnv_in_known_regions_correct.tsv",
        sep="\t",
    )

    df_combined_nt_correct = pd.read_csv(
        "tests/test_files/df_combined_nt_correct.tsv", sep="\t"
    )

    def test_find_exact_matches(self):
        df_exact_matches = find_exact_matches(
            df_mutations=df_mutations,
            df_resistance_variants=df_nt_resistance_variants,
            dict_col_rename=dict_col_rename_nt,
        )
        self.assertEqual(df_exact_matches.shape[0], 1)
        self.assertEqual(df_exact_matches.shape[1], 11)
        pd.testing.assert_frame_equal(df_exact_matches, self.df_exact_matches_correct)

    def test_find_large_indels_should_find_novel(self):
        df_large_indels = find_large_indels(
            df_mutations=df_mutations,
            screen_region=(350, 450),
            genetic_element="b0004_promoter",
        )
        self.assertEqual(df_large_indels.shape[0], 1)
        self.assertEqual(df_large_indels.shape[1], 10)
        df_large_indels.reset_index(drop=True, inplace=True)
        self.df_find_large_indels_correct.reset_index(drop=True, inplace=True)
        self.assertTrue(df_large_indels.equals(self.df_find_large_indels_correct))

    def test_find_large_indels_should_not_find(self):
        df_large_indels = find_large_indels(
            df_mutations=df_mutations,
            screen_region=(450, 550),
            genetic_element="b0005_promoter",
        )
        # Just check if empty
        self.assertEqual(df_large_indels.shape[0], 0)
        self.assertEqual(df_large_indels.shape[1], 10)

    def test_screen_for_possible_cnv_in_known_regions(self):
        df_screen_for_possible_cnv_in_known_regions = (
            screen_for_possible_cnv_in_known_regions(
                df_mutations=df_mutations,
                df_resistance_variants=df_nt_resistance_variants,
                dict_col_rename=dict_col_rename_nt,
            )
        )
        # find two hits (one should be filtered out later)
        self.assertEqual(df_screen_for_possible_cnv_in_known_regions.shape[0], 2)
        self.assertEqual(df_screen_for_possible_cnv_in_known_regions.shape[1], 11)
        df_screen_for_possible_cnv_in_known_regions.reset_index(drop=True, inplace=True)
        self.df_screen_for_possible_cnv_in_known_regions_correct.reset_index(
            drop=True, inplace=True
        )
        self.assertTrue(
            df_screen_for_possible_cnv_in_known_regions.equals(
                self.df_screen_for_possible_cnv_in_known_regions_correct
            )
        )

    def test_combine_exact_matches_and_possible_cnvs(self):
        # Copy to not change in memory dfs
        df_combined = combine_exact_matches_and_possible_cnvs(
            df_exact_matches=self.df_exact_matches_correct.copy(),
            df_possible_cnvs=self.df_screen_for_possible_cnv_in_known_regions_correct.copy(),
        )
        df_combined.to_csv(
            "tests/test_files/df_combined_nt_correct.tsv", sep="\t", index=False
        )

        self.assertEqual(df_combined.shape[0], 2)
        self.assertEqual(df_combined.shape[1], 11)

        df_combined.reset_index(drop=True, inplace=True)
        self.df_exact_matches_correct.reset_index(drop=True, inplace=True)
        self.assertTrue(df_combined.equals(self.df_combined_nt_correct))

