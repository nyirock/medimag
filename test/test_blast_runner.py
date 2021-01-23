import os
import traceback
import pandas as pd
from unittest import TestCase

from pandas.testing import assert_frame_equal

from tools.blast_runner import BlastRunner
from tools.fs import init_blast_dir

TEST_DATA_DIR = 'data'
BLAST_DIR_NAME = 'blast'
BLAST_DB_DIR_NAME = "blast_db/db"

#TODO check for newer blast version and test accordingly
class TestBlastRunner(TestCase):
    pmoA_reference_path = os.path.join(TEST_DATA_DIR, "all_pmoA_nr99.fasta")
    db_path = None
    def __init__(self, *args, **kwargs):
        super(TestBlastRunner, self).__init__(*args, **kwargs)


        #self.run_op_finders(self)
    @classmethod
    def setUpClass(cls):
        cls.db_path = init_blast_dir(TEST_DATA_DIR, BLAST_DIR_NAME, BLAST_DB_DIR_NAME)
        # cls.init_out_fnames(cls.input_basename, cls.out, cls.testdata_dir)
        # cls.run_op_finder(cls.args, cls.python_path)



    def setUp(self) -> None:
        self.large_mg_fasta_path = os.path.join(TEST_DATA_DIR, 'mg_large.fasta')
        self.large_mg_blast_out_path = os.path.join(TEST_DATA_DIR, BLAST_DIR_NAME, 'mg_large.tab')
        self.large_mg_blast_expected_name = os.path.join(TEST_DATA_DIR, "recruited_mg_0.tab")

        self.no_pmoA_hits_fasta_path = os.path.join(TEST_DATA_DIR, "mg_no_pmoA_hits.fasta")
        self.no_pmoA_blast_out_path = os.path.join(TEST_DATA_DIR, "mg_no_pmoA_hits.tab")

        try:
            self.mg_large_blast_expected = pd.read_csv(self.large_mg_blast_expected_name, sep='\t', header=None).sort_values(by=[0, 1], ignore_index=True)
        except IOError:
            tb = traceback.format_exc()
            print(f"Error while loading the fixture file.")
            print(tb)



    def test_run_blast_parallel_large_non_empty(self):
        blast_runner = BlastRunner(self.__class__.pmoA_reference_path, self.__class__.db_path)
        blast_runner.run_blast_parallel(self.large_mg_fasta_path, self.large_mg_blast_out_path)
        mg_large_blast_df = pd.read_csv(self.large_mg_blast_out_path, sep='\t', header=None).sort_values(by=[0, 1], ignore_index=True)
        assert_frame_equal(self.mg_large_blast_expected, mg_large_blast_df)

    def test_run_blast_parallel_large_empty(self):
        blast_runner = BlastRunner(self.__class__.pmoA_reference_path, self.__class__.db_path)
        blast_runner.run_blast_parallel(self.no_pmoA_hits_fasta_path, self.no_pmoA_blast_out_path)
        with self.assertRaises(pd.errors.EmptyDataError):
            pd.read_csv(self.no_pmoA_blast_out_path)