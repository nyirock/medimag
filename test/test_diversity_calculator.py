import os
import traceback
from unittest import TestCase
from pandas.testing import assert_frame_equal
import pandas as pd

from tools.diversity_calculator import DiversityCalculator

TEST_DATA_DIR = 'data'


class TestDiversityCalculator(TestCase):

    def setUp(self) -> None:
        self.single_file_blast_runner_csv0 = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.csv')
        self.single_file_blast_runner_csv1 = os.path.join(TEST_DATA_DIR, 'recruited_mg_1.csv')
        self.single_file_blast_runner_csv2 = os.path.join(TEST_DATA_DIR, 'recruited_mg_2.csv')
        self.m7_1y2_r12unpaired_csv3 = os.path.join(TEST_DATA_DIR, 'm7_1y2_r12unpaired.csv')

        self.single_file_blast_runner_csv_empty = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.csv')
        self.sigle_file_blast_runner_csv_loaded_w_frac = os.path.join(TEST_DATA_DIR, 'recruited_mg_0_w_frac.pkl')
        self.sigle_file_blast_runner_csv_loaded = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.pkl')
        self.three_df_blast_runner_loaded_w_frac = os.path.join(TEST_DATA_DIR, 'merged_3_dfs_w_frac.pkl')

        try:
            self.single_file_data0 = pd.read_csv(self.single_file_blast_runner_csv0, sep='\t').drop("Unnamed: 0",
                                                                                                    axis=1)
            self.data0_hit_summary = pd.read_pickle(os.path.join(TEST_DATA_DIR, "m7_1y2_r1_hit_summary.pkl"))
            self.data_m7_1y2_r12unpaired_hit_summary = pd.read_pickle(os.path.join(TEST_DATA_DIR, "m7_1y2_r12unpaired_hit_summary.pkl"))

            self.single_file_data1 = pd.read_csv(self.single_file_blast_runner_csv1, sep='\t').drop("Unnamed: 0",
                                                                                                    axis=1)
            self.single_file_data2 = pd.read_csv(self.single_file_blast_runner_csv2, sep='\t').drop("Unnamed: 0",
                                                                                                    axis=1)
            # self.single_file_data_empty = pd.read_csv(self.single_file_blast_runner_csv_empty, sep='\t')
            self.sinlge_file_load_result_w_frac = pd.read_pickle(self.sigle_file_blast_runner_csv_loaded_w_frac).drop(
                "Unnamed: 0", axis=1)
            self.sinlge_file_load_result = pd.read_pickle(self.sigle_file_blast_runner_csv_loaded).drop("Unnamed: 0",
                                                                                                        axis=1)
            self.three_df_file_load_result_w_frac = pd.read_pickle(self.three_df_blast_runner_loaded_w_frac).drop(
                "Unnamed: 0", axis=1)

        except IOError:
            tb = traceback.format_exc()
            print(f"Error while loading the fixture file.")
            print(tb)

        self.dc0 = DiversityCalculator(("", self.single_file_blast_runner_csv0), TEST_DATA_DIR)
        self.dc_m7_1y2 = DiversityCalculator(("", self.m7_1y2_r12unpaired_csv3), TEST_DATA_DIR)

    ####### __init Tests #############
    def test___init_no_blast_hits(self):
        sname = "name"
        dc = DiversityCalculator((sname, self.single_file_blast_runner_csv_empty), TEST_DATA_DIR)
        self.assertTrue(dc.reference_hits.empty)
        self.assertTrue(dc.no_hits)
        self.assertEqual(sname, dc.sname)

    def test___init_single_non_empty(self):
        sname = "name"
        dc = DiversityCalculator((sname, self.single_file_blast_runner_csv0), TEST_DATA_DIR)
        assert_frame_equal(dc.reference_hits, self.sinlge_file_load_result_w_frac)
        self.assertFalse(dc.reference_hits.empty)
        self.assertFalse(dc.no_hits)
        self.assertEqual(sname, dc.sname)

    # def test___init_three_files_plus_empty(self):
    #     sname = "name"
    #     dc = DiversityCalculator([self.single_file_blast_runner_csv0,
    #                               self.single_file_blast_runner_csv1,
    #                               self.single_file_blast_runner_csv2,
    #                               self.single_file_blast_runner_csv_empty], sname)
    #     assert_frame_equal(dc.reference_hits, self.three_df_file_load_result_w_frac)
    #     self.assertFalse(dc.reference_hits.empty)
    #     self.assertFalse(dc.no_hits)
    #     self.assertEqual(sname, dc.sample_name)

    ####### _load_hit_summary Tests #############
    def test__load_hits_single_file(self):
        raw_hit_summary = DiversityCalculator._load_hits(self.single_file_blast_runner_csv0)
        raw_hit_summary['iden'] = raw_hit_summary['iden'] / 100
        assert_frame_equal(raw_hit_summary, self.sinlge_file_load_result_w_frac)


    def test__load_hits_empty_file(self):
        raw_hit_summary = DiversityCalculator._load_hits(self.single_file_blast_runner_csv_empty)
        self.assertTrue(raw_hit_summary.empty)

    # def test__load_hits_3_files(self):
    #     raw_hit_summary = DiversityCalculator._load_hits([self.single_file_blast_runner_csv0,
    #                                                       self.single_file_blast_runner_csv1,
    #                                                       self.single_file_blast_runner_csv2, ])
    #     # with self.assertRaises(AssertionError):
    #     #     assert_frame_equal(raw_hit_summary, self.sinlge_file_load_result_w_frac)
    #     raw_hit_summary.iden = raw_hit_summary.iden / 100
    #     assert_frame_equal(raw_hit_summary, self.three_df_file_load_result_w_frac)

    # def test__load_hits_3_files_and_empty(self):
    #     raw_hit_summary1 = DiversityCalculator._load_hits([self.single_file_blast_runner_csv0,
    #                                                        self.single_file_blast_runner_csv1,
    #                                                        self.single_file_blast_runner_csv2,
    #                                                        self.single_file_blast_runner_csv_empty])
    #     # with self.assertRaises(AssertionError):
    #     #     assert_frame_equal(raw_hit_summary, self.sinlge_file_load_result_w_frac)
    #     raw_hit_summary1.iden = raw_hit_summary1.iden / 100
    #     #self.assertTrue((raw_hit_summary1 == self.three_df_file_load_result_w_frac).all().all())
    #     assert_frame_equal(raw_hit_summary1, self.three_df_file_load_result_w_frac, check_dtype=False)
    #     assert_frame_equal(raw_hit_summary1, self.three_df_file_load_result_w_frac)

    ####### _convert_iden_percent_to_frac Tests #############
    def test__convert_iden_percent_to_frac_single_file(self):

        assert_frame_equal(self.dc0.reference_hits, self.sinlge_file_load_result_w_frac)

    # def test__convert_iden_percent_to_frac_three_files(self):
    #     dc = DiversityCalculator([self.single_file_blast_runner_csv0,
    #                               self.single_file_blast_runner_csv1,
    #                               self.single_file_blast_runner_csv2, ], "")
    #     assert_frame_equal(dc.reference_hits, self.three_df_file_load_result_w_frac)

    def test__summarize_reference_hits(self):
        assert_frame_equal(self.data0_hit_summary, self.dc0.reference_hit_summary)

    def test__summarize_reference_hits2(self):
        assert_frame_equal(self.data_m7_1y2_r12unpaired_hit_summary, self.dc_m7_1y2.reference_hit_summary)

    def test__calculate_concordance_stats(self):
        exp_conc_stats = (0.9599823028215652, 0.07893666529125956, 0.013537526413318799)
        for exp, actual in zip(exp_conc_stats, (self.dc0.concordance.mean, self.dc0.concordance.std, self.dc0.concordance.sem)):
            with self.subTest(msg=f"{exp}: {actual}"):
                self.assertAlmostEqual(exp, actual)

    def test__calculate_concordance_stats2(self):
        exp_conc_stats = (0.9648019939042705, 0.06604612828304085, 0.00800926950217857)
        for exp, actual in zip(exp_conc_stats, (self.dc_m7_1y2.concordance.mean, self.dc_m7_1y2.concordance.std, self.dc_m7_1y2.concordance.sem)):
            with self.subTest(msg=f"{exp}: {actual}"):
                self.assertAlmostEqual(exp, actual)

