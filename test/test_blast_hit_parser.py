import os
import traceback
from unittest import TestCase
import pandas as pd
from pandas._testing import assert_frame_equal


from tools.blast_hit_parser import BlastHitParser
from tools.fasta_utils import parse_contigs_ind
from tools.fs import load_blast_hits_to_df

TEST_DATA_DIR = 'data'
REFERENCE = os.path.join(TEST_DATA_DIR, "all_pmoA_nr99.fasta")

blast_out_paths = [os.path.join(TEST_DATA_DIR, 'recruited_mg_0.tab'),
                   os.path.join(TEST_DATA_DIR, 'recruited_mg_0_1.tab')]

class TestBlastHitParser(TestCase):
    ref_index = None
    def __init__(self, *args, **kwargs):
        super(TestBlastHitParser, self).__init__(*args, **kwargs)

    @classmethod
    def setUpClass(cls):
        cls.ref_index = parse_contigs_ind(REFERENCE)


    def setUp(self) -> None:
        #self.pmoA_reference_path = os.path.join(TEST_DATA_DIR, "all_pmoA_nr99.fasta")
        self.mg_with_pmoA_hits_fasta_path = os.path.join(TEST_DATA_DIR, "mg_large.fasta")
        self.mg_no_pmoA_hits_fasta_path = os.path.join(TEST_DATA_DIR, "mg_no_pmoA_hits.fasta")

        # self.blast_out_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.tab')
        # self.blast_out_path2 = os.path.join(TEST_DATA_DIR, 'recruited_mg_0_1.tab')
        self.blast_out_parsed_expected_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.csv')
        self.blast_out_parsed_expected_fasta_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.fasta')

        self.blast_out_no_hits_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.tab')
        self.blast_out_parsed_no_hits_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.csv')
        self.blast_out_parsed_no_hits_expected_fasta_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.fasta')

        try:
            self.blast_out_expected = pd.read_csv(self.blast_out_parsed_expected_path, sep='\t').sort_values(by=['quid', 'suid'], ignore_index=True).drop("Unnamed: 0", axis=1)
            self.blast_out_no_hits_expected = pd.read_csv(self.blast_out_parsed_no_hits_path, sep='\t').sort_values(by=['quid', 'suid'], ignore_index=True).drop("Unnamed: 0",axis=1)
        except IOError:
            tb = traceback.format_exc()
            print(f"Error while loading the fixture file.")
            print(tb)


    def test_parse_blast_hits_not_empty(self):
        for blast_out_path in blast_out_paths:
            with self.subTest(msg=blast_out_path):
                sample_name = "test_blast_parser1"
                blast_df = load_blast_hits_to_df(blast_out_path)
                bhp = BlastHitParser((sample_name, self.mg_with_pmoA_hits_fasta_path ),
                                     self.__class__.ref_index, blast_df, TEST_DATA_DIR)
                bhp.parse_blast_hits()
                parsed_hits_df = pd.read_csv(os.path.join(TEST_DATA_DIR, sample_name+".csv"),
                                             sep='\t').drop("Unnamed: 0",axis=1).sort_values(by=['quid', 'suid'], ignore_index=True)
                assert_frame_equal(self.blast_out_expected, parsed_hits_df)
                #TODO: make a better way of comparing fasta files than by size
                self.assertEqual(os.path.getsize(os.path.join(TEST_DATA_DIR, sample_name+".fasta")),
                                 os.path.getsize(self.blast_out_parsed_expected_fasta_path))

    def test_parse_blast_hits_empty(self):

        sample_name = "test_blast_parser2"
        blast_df = load_blast_hits_to_df(self.blast_out_no_hits_path)
        bhp = BlastHitParser((sample_name, self.mg_no_pmoA_hits_fasta_path),
                             self.__class__.ref_index, blast_df, TEST_DATA_DIR)
        bhp.parse_blast_hits()
        parsed_hits_df = pd.read_csv(os.path.join(TEST_DATA_DIR, sample_name + ".csv"),
                                     sep='\t').drop("Unnamed: 0", axis=1).sort_values(by=['quid', 'suid'],
                                                                                      ignore_index=True)
        assert_frame_equal(self.blast_out_no_hits_expected, parsed_hits_df)
        # TODO: make a better way of comparing fasta files than by size
        self.assertEqual(os.path.getsize(os.path.join(TEST_DATA_DIR, sample_name + ".fasta")),
                         os.path.getsize(self.blast_out_parsed_no_hits_expected_fasta_path))