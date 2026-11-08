"""Tests for the parsing functions of the G2P, Collins et al. 2022, GeneReviews, Open Targets gene burden, and
ACMG secondary findings downloads, and for the HGNC gene symbol lookup they rely on.

Run from the repo root with:  python3 -m unittest annotation_utils/gene_annotation_sources_tests.py
"""

import gzip
import os
import tempfile
import unittest

import pandas as pd

import annotation_utils.cache_utils as cache_utils
from annotation_utils.get_acmg_secondary_findings_genes import (
    convert_acmg_secondary_findings_response_to_table, parse_acmg_secondary_findings_list_version)
from annotation_utils.get_collins_dosage_sensitivity_scores import parse_collins_dosage_sensitivity_scores
from annotation_utils.get_g2p_table import convert_g2p_records_to_gene_disease_table, find_latest_g2p_release_date
from annotation_utils.get_genereviews_table import combine_genereviews_gene_symbols_and_titles
from annotation_utils.get_hgnc_table import build_gene_symbol_to_ensg_id_map
from annotation_utils.get_open_targets_gene_burden_table import (
    convert_gene_burden_evidence_to_table, find_latest_open_targets_release, find_parquet_file_names)


class ReturnStaleCachedTableOrRaiseTests(unittest.TestCase):

    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.original_cache_dir = cache_utils.CACHE_DIR
        cache_utils.CACHE_DIR = self.temporary_directory.name

    def tearDown(self):
        cache_utils.CACHE_DIR = self.original_cache_dir
        self.temporary_directory.cleanup()

    def test_returns_the_cached_table_when_one_exists(self):
        pd.DataFrame({"gene_id": ["ENSG1"]}).to_csv(
            os.path.join(self.temporary_directory.name, "g2p_table.abc123.tsv.gz"), sep="\t", index=False)
        df = cache_utils.return_stale_cached_table_or_raise("get_g2p_table", "G2P", ValueError("download failed"))
        self.assertEqual(df["gene_id"].tolist(), ["ENSG1"])
        self.assertTrue(df.attrs[cache_utils.FROM_CACHE_ATTR])

    def test_reraises_the_download_error_when_nothing_is_cached(self):
        with self.assertRaisesRegex(ValueError, "download failed"):
            cache_utils.return_stale_cached_table_or_raise("get_g2p_table", "G2P", ValueError("download failed"))


class BuildGeneSymbolToEnsgIdMapTests(unittest.TestCase):

    def setUp(self):
        self.df_hgnc = pd.DataFrame({
            "Approved symbol": ["GENE1", "GENE2", "GENE3", "GENE4", "NOENSG"],
            "Previous symbols": ["OLD1, SHARED_OLD", "GENE1_LIKE, SHARED_OLD", "GENE2", None, "OLD_NOENSG"],
            "Alias symbols": ["ALIAS1", "OLD1", "ALIAS3, ALIAS_DUP", "ALIAS_DUP", "ALIAS_NOENSG"],
            "Ensembl gene ID": ["ENSG1", "ENSG2", "ENSG3", "ENSG4", None],
        })

    def test_approved_symbol_wins_over_another_genes_previous_symbol(self):
        # GENE2 is GENE3's previous symbol, but it is GENE2's approved symbol
        self.assertEqual(build_gene_symbol_to_ensg_id_map(self.df_hgnc)["GENE2"], "ENSG2")

    def test_unique_previous_symbol_and_unique_alias_resolve(self):
        symbol_to_ensg_id = build_gene_symbol_to_ensg_id_map(self.df_hgnc)
        self.assertEqual(symbol_to_ensg_id["GENE1_LIKE"], "ENSG2")
        self.assertEqual(symbol_to_ensg_id["ALIAS3"], "ENSG3")

    def test_previous_symbol_wins_over_another_genes_alias(self):
        # OLD1 is GENE1's previous symbol and GENE2's alias
        self.assertEqual(build_gene_symbol_to_ensg_id_map(self.df_hgnc)["OLD1"], "ENSG1")

    def test_ambiguous_symbols_and_genes_without_ensg_ids_are_left_out(self):
        symbol_to_ensg_id = build_gene_symbol_to_ensg_id_map(self.df_hgnc)
        for symbol in ["SHARED_OLD", "ALIAS_DUP", "NOENSG", "OLD_NOENSG", "ALIAS_NOENSG"]:
            self.assertNotIn(symbol, symbol_to_ensg_id)


class G2PTests(unittest.TestCase):

    def test_find_latest_g2p_release_date(self):
        listing = '<a href="2026_09_28/">2026_09_28/</a> <a href="2026_10_28/">x</a> <a href="Data_download_format_202511.txt">'
        self.assertEqual(find_latest_g2p_release_date(listing), "2026_10_28")
        with self.assertRaises(ValueError):
            find_latest_g2p_release_date("<html></html>")

    def test_convert_g2p_records_to_gene_disease_table(self):
        df_g2p = pd.DataFrame({
            "g2p id": ["G2P00001", "G2P00002"],
            "hgnc id": [5017, 663],
            "disease name": ["HMX1-related oculoauricular syndrome", "ABC-related disorder"],
            "allelic requirement": ["biallelic_autosomal", "some_new_requirement"],
            "confidence": ["definitive", None],
            "molecular mechanism": ["loss of function", "gain of function"],
            "variant consequence": ["absent gene product", "absent gene product; altered gene product structure"],
            "panel": ["DD; Eye", "Skin"],
        })
        df = convert_g2p_records_to_gene_disease_table(df_g2p)
        self.assertEqual(df["hgnc_gene_id"].tolist(), ["HGNC:5017", "HGNC:663"])
        self.assertEqual(df["inheritance"].tolist(), ["AR", "some_new_requirement"])
        self.assertEqual(df["confidence"].tolist()[0], "definitive")
        self.assertTrue(pd.isna(df["confidence"].tolist()[1]))
        self.assertEqual(df["variant_consequence"].tolist()[1], "absent gene product, altered gene product structure")
        self.assertEqual(df["panel"].tolist()[0], "DD, Eye")


class CollinsDosageSensitivityScoresTests(unittest.TestCase):

    def test_parse_collins_dosage_sensitivity_scores(self):
        file_contents = gzip.compress(b"#gene\tpHaplo\tpTriplo\nCACNA1C\t0.998\t1\nZNF462\t1\t0.988\n")
        df = parse_collins_dosage_sensitivity_scores(file_contents)
        self.assertEqual(df.columns.tolist(), ["gene_symbol", "pHaplo", "pTriplo"])
        self.assertEqual(df["gene_symbol"].tolist(), ["CACNA1C", "ZNF462"])
        self.assertAlmostEqual(df["pTriplo"].tolist()[1], 0.988)

    def test_duplicate_gene_symbols_raise(self):
        file_contents = gzip.compress(b"#gene\tpHaplo\tpTriplo\nCACNA1C\t0.998\t1\nCACNA1C\t1\t0.988\n")
        with self.assertRaises(ValueError):
            parse_collins_dosage_sensitivity_scores(file_contents)


class GeneReviewsTests(unittest.TestCase):

    def setUp(self):
        self.df_titles = pd.DataFrame({
            "#GR_shortname": ["trimethylaminuria", "cdls"],
            "GR_Title": ["Primary Trimethylaminuria", "Cornelia de Lange Syndrome"],
            "NBK_id": ["NBK1103", "NBK1104"],
            "PMID": [20301282, 0],
        })

    def test_combine_genereviews_gene_symbols_and_titles(self):
        df_gene_symbols = pd.DataFrame({
            "#NBK_id": ["NBK1103", None, "NBK1104", "NBK1104", "NBK1104"],
            "GR_shortname": ["trimethylaminuria", "cln2-dis", "cdls", "cdls", "cdls"],
            "genesymbol": ["FMO3", "TPP1", "BRD4", "HDAC8", "Not applicable"],
        })
        df = combine_genereviews_gene_symbols_and_titles(df_gene_symbols, self.df_titles)
        self.assertEqual(df.to_dict("records"), [
            {"gene_symbol": "FMO3", "nbk_id": "NBK1103", "title": "Primary Trimethylaminuria"},
            {"gene_symbol": "BRD4", "nbk_id": "NBK1104", "title": "Cornelia de Lange Syndrome"},
            {"gene_symbol": "HDAC8", "nbk_id": "NBK1104", "title": "Cornelia de Lange Syndrome"},
        ])

    def test_chapter_without_title_raises(self):
        df_gene_symbols = pd.DataFrame({"#NBK_id": ["NBK9999"], "GR_shortname": ["x"], "genesymbol": ["GENE1"]})
        with self.assertRaises(ValueError):
            combine_genereviews_gene_symbols_and_titles(df_gene_symbols, self.df_titles)


class OpenTargetsGeneBurdenTests(unittest.TestCase):

    def test_find_latest_open_targets_release(self):
        listing = '<a href="9.12/"> <a href="25.12/"> <a href="26.09/"> <a href="26.03/"> <a href="latest/"> <a href="master/">'
        self.assertEqual(find_latest_open_targets_release(listing), "26.09")
        with self.assertRaises(ValueError):
            find_latest_open_targets_release('<a href="latest/">')

    def test_find_parquet_file_names(self):
        listing = '<a href="/pub/">..</a> <a href="_SUCCESS"> <a href="part-00001-b.zstd.parquet"> <a href="part-00000-a.zstd.parquet">'
        self.assertEqual(find_parquet_file_names(listing), ["part-00000-a.zstd.parquet", "part-00001-b.zstd.parquet"])

    def test_convert_gene_burden_evidence_to_table(self):
        df_evidence = pd.DataFrame({
            "targetId": ["ENSG00000118137", "ENSG2", "ENSG3", "ENSG4"],
            "diseaseId": ["EFO_0022282", "MONDO_0005258", "MONDO_0005258", "MONDO_0005258"],
            "diseaseFromSource": ["HDL trait", "Autism", "Autism", "Autism"],
            "projectId": ["AstraZeneca PheWAS Portal", None, None, None],
            "cohortId": ["UK Biobank 470k", "UK Biobank, MDCS,  MCPS", None, None],
            "literature": [["34375979"], None, ["35982159"], None],
            "statisticalMethod": ["flexdmg", "TADA", "TADA", "TADA"],
            "pValueMantissa": [1.677, 2.0, 1.0, 1.0],
            "pValueExponent": [-15, -400, -6, -6],
            "directionOnTrait": ["protect", None, "risk", "risk"],
        })
        df = convert_gene_burden_evidence_to_table(df_evidence)
        self.assertEqual(df["disease_id"].tolist(), ["EFO:0022282", "MONDO:0005258", "MONDO:0005258", "MONDO:0005258"])
        self.assertEqual(df["study"].tolist(), ["AstraZeneca PheWAS Portal", "UK Biobank, MDCS, MCPS", "PMID:35982159", ""])
        self.assertAlmostEqual(df["p_value"].tolist()[0], 1.677e-15, delta=1e-27)
        self.assertEqual(df["p_value"].tolist()[1], 0.0)
        self.assertEqual(df["direction_on_trait"].tolist(), ["protect", None, "risk", "risk"])
        self.assertEqual(df.columns.tolist(), ["gene_id", "disease_id", "study", "p_value", "direction_on_trait"])


class AcmgSecondaryFindingsTests(unittest.TestCase):

    def test_parse_acmg_secondary_findings_list_version(self):
        self.assertEqual(parse_acmg_secondary_findings_list_version("<h3>ACMG SF v3.3 Curations</h3> (SF v3.2)"), "3.3")
        self.assertEqual(parse_acmg_secondary_findings_list_version("<html>no version</html>"), "")

    def test_convert_acmg_secondary_findings_response_to_table(self):
        response_json = {"total": 2, "rows": [{"symbol": "ABCD1", "hgnc_id": "HGNC:61"}, {"symbol": "PLN", "hgnc_id": "HGNC:9080"}]}
        df = convert_acmg_secondary_findings_response_to_table(response_json)
        self.assertEqual(df.to_dict("records"), [
            {"hgnc_gene_id": "HGNC:61", "gene_symbol": "ABCD1"},
            {"hgnc_gene_id": "HGNC:9080", "gene_symbol": "PLN"},
        ])

    def test_incomplete_response_raises(self):
        with self.assertRaises(ValueError):
            convert_acmg_secondary_findings_response_to_table({"total": 84, "rows": [{"symbol": "ABCD1", "hgnc_id": "HGNC:61"}]})


if __name__ == "__main__":
    unittest.main()
