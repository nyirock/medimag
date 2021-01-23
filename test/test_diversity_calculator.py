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
        self.single_file_blast_runner_csv_empty = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.csv')
        self.sigle_file_blast_runner_csv_loaded_w_frac = os.path.join(TEST_DATA_DIR, 'recruited_mg_0_w_frac.pkl')
        self.sigle_file_blast_runner_csv_loaded = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.pkl')
        self.three_df_blast_runner_loaded_w_frac = os.path.join(TEST_DATA_DIR, 'merged_3_dfs_w_frac.pkl')
        try:
            self.single_file_data0 = pd.read_csv(self.single_file_blast_runner_csv0, sep='\t').drop("Unnamed: 0",
                                                                                                    axis=1)
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

    ####### __init Tests #############
    def test___init_no_blast_hits(self):
        sname = "name"
        dc = DiversityCalculator([self.single_file_blast_runner_csv_empty], sname)
        self.assertTrue(dc.reference_hits.empty)
        self.assertTrue(dc.no_hits)
        self.assertEqual(sname, dc.sample_name)

    def test___init_single_non_empty(self):
        sname = "name"
        dc = DiversityCalculator([self.single_file_blast_runner_csv0], sname)
        assert_frame_equal(dc.reference_hits, self.sinlge_file_load_result_w_frac)
        self.assertFalse(dc.reference_hits.empty)
        self.assertFalse(dc.no_hits)
        self.assertEqual(sname, dc.sample_name)

    def test___init_three_files_plus_empty(self):
        sname = "name"
        dc = DiversityCalculator([self.single_file_blast_runner_csv0,
                                  self.single_file_blast_runner_csv1,
                                  self.single_file_blast_runner_csv2,
                                  self.single_file_blast_runner_csv_empty], sname)
        assert_frame_equal(dc.reference_hits, self.three_df_file_load_result_w_frac)
        self.assertFalse(dc.reference_hits.empty)
        self.assertFalse(dc.no_hits)
        self.assertEqual(sname, dc.sample_name)

    ####### _load_hit_summary Tests #############
    def test__load_hits_single_file(self):
        raw_hit_summary = DiversityCalculator._load_hits([self.single_file_blast_runner_csv0])

        with self.assertRaises(AssertionError):
            assert_frame_equal(raw_hit_summary, self.sinlge_file_load_result_w_frac)

        assert_frame_equal(raw_hit_summary, self.sinlge_file_load_result)

    def test__load_hits_empty_file(self):
        raw_hit_summary = DiversityCalculator._load_hits([self.single_file_blast_runner_csv_empty])
        self.assertTrue(raw_hit_summary.empty)

    def test__load_hits_3_files(self):
        raw_hit_summary = DiversityCalculator._load_hits([self.single_file_blast_runner_csv0,
                                                          self.single_file_blast_runner_csv1,
                                                          self.single_file_blast_runner_csv2, ])
        # with self.assertRaises(AssertionError):
        #     assert_frame_equal(raw_hit_summary, self.sinlge_file_load_result_w_frac)
        raw_hit_summary.iden = raw_hit_summary.iden / 100
        assert_frame_equal(raw_hit_summary, self.three_df_file_load_result_w_frac)

    def test__load_hits_3_files_and_empty(self):
        raw_hit_summary1 = DiversityCalculator._load_hits([self.single_file_blast_runner_csv0,
                                                           self.single_file_blast_runner_csv1,
                                                           self.single_file_blast_runner_csv2,
                                                           self.single_file_blast_runner_csv_empty])
        # with self.assertRaises(AssertionError):
        #     assert_frame_equal(raw_hit_summary, self.sinlge_file_load_result_w_frac)
        raw_hit_summary1.iden = raw_hit_summary1.iden / 100
        #self.assertTrue((raw_hit_summary1 == self.three_df_file_load_result_w_frac).all().all())
        assert_frame_equal(raw_hit_summary1, self.three_df_file_load_result_w_frac, check_dtype=False)
        assert_frame_equal(raw_hit_summary1, self.three_df_file_load_result_w_frac)

    ####### _convert_iden_percent_to_frac Tests #############
    def test__convert_iden_percent_to_frac_single_file(self):
        dc = DiversityCalculator([self.single_file_blast_runner_csv0], "")
        assert_frame_equal(dc.reference_hits, self.sinlge_file_load_result_w_frac)

    def test__convert_iden_percent_to_frac_three_files(self):
        dc = DiversityCalculator([self.single_file_blast_runner_csv0,
                                  self.single_file_blast_runner_csv1,
                                  self.single_file_blast_runner_csv2, ], "")
        assert_frame_equal(dc.reference_hits, self.three_df_file_load_result_w_frac)

