import os
import traceback
from unittest import TestCase
from pandas._testing import assert_frame_equal
import pandas as pd

from tools.blast_hit_parser_fast import BlastHitParserFast
from tools.fasta_utils import parse_contigs_ind

TEST_DATA_DIR = 'data'
REFERENCE = os.path.join(TEST_DATA_DIR, "all_pmoA_nr99.fasta")

blast_out_paths = [os.path.join(TEST_DATA_DIR, 'recruited_mg_0_w_qlen.tab'),
                   ]


class TestBlastHitParserFast(TestCase):
    ref_index = None

    def __init__(self, *args, **kwargs):
        super(TestBlastHitParserFast, self).__init__(*args, **kwargs)

    @classmethod
    def setUpClass(cls):
        cls.ref_index = parse_contigs_ind(REFERENCE)

    def setUp(self) -> None:
        self.mg_with_pmoA_hits_fasta_path = os.path.join(TEST_DATA_DIR, "mg_large.fasta")
        self.mg_no_pmoA_hits_fasta_path = os.path.join(TEST_DATA_DIR, "mg_no_pmoA_hits.fasta")

        self.mg_with_pmoA_hits_fasta_path_gz = os.path.join(TEST_DATA_DIR, "mg_large.fasta.gz")
        self.mg_no_pmoA_hits_fasta_path_gz = os.path.join(TEST_DATA_DIR, "mg_no_pmoA_hits.fasta.gz")

        self.blast_out_parsed_expected_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.csv')
        self.blast_out_parsed_expected_fasta_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_0.fasta')

        self.blast_out_no_hits_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.tab')
        self.blast_out_parsed_no_hits_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.csv')
        self.blast_out_parsed_no_hits_expected_fasta_path = os.path.join(TEST_DATA_DIR, 'recruited_mg_8.fasta')

        try:
            self.blast_out_parsed_expected = pd.read_csv(self.blast_out_parsed_expected_path, sep='\t').sort_values(by=['quid', 'suid'], ignore_index=True).drop("Unnamed: 0", axis=1)
            self.blast_out_no_hits_expected = pd.read_csv(self.blast_out_parsed_no_hits_path, sep='\t').sort_values(by=['quid', 'suid'], ignore_index=True).drop("Unnamed: 0",axis=1)
        except IOError:
            tb = traceback.format_exc()
            print(f"Error while loading the fixture file.")
            print(tb)

    def test_parse_blast_hits_not_empty(self):
        for blast_out_path in blast_out_paths:
            with self.subTest(msg=blast_out_path):
                sample_name = "test_blast_parser_fast1"
                #blast_df = load_blast_hits_to_df(blast_out_path)
                bhp = BlastHitParserFast((sample_name, self.mg_with_pmoA_hits_fasta_path ),
                                     self.__class__.ref_index, blast_out_path, TEST_DATA_DIR, outdir_name="")
                out_path = bhp.parse_blast_hits()
                expected_df = self.blast_out_parsed_expected.drop(['Contig_nt', 'Contig_GC'], axis=1)
                parsed_hits_df = pd.read_csv(out_path,
                                             sep='\t').drop("Unnamed: 0",axis=1).sort_values(by=['quid', 'suid'], ignore_index=True)
                assert_frame_equal(expected_df, parsed_hits_df)

    def test_parse_blast_hits_not_empty_gz(self):
        for blast_out_path in blast_out_paths:
            with self.subTest(msg=blast_out_path):
                sample_name = "test_blast_parser_fast1"
                #blast_df = load_blast_hits_to_df(blast_out_path)
                bhp = BlastHitParserFast((sample_name, self.mg_with_pmoA_hits_fasta_path_gz),
                                     self.__class__.ref_index, blast_out_path, TEST_DATA_DIR, outdir_name="")
                out_path = bhp.parse_blast_hits()
                expected_df = self.blast_out_parsed_expected.drop(['Contig_nt', 'Contig_GC'], axis=1)
                parsed_hits_df = pd.read_csv(out_path,
                                             sep='\t').drop("Unnamed: 0",axis=1).sort_values(by=['quid', 'suid'], ignore_index=True)
                assert_frame_equal(expected_df, parsed_hits_df)

    def test_parse_blast_hits_empty(self):
        file_lst = [self.mg_no_pmoA_hits_fasta_path, self.mg_no_pmoA_hits_fasta_path_gz]
        for file_path in file_lst:
            with self.subTest(msg=file_path):
                sample_name = "test_blast_parser2"
                #blast_df = load_blast_hits_to_df(self.blast_out_no_hits_path)
                bhp = BlastHitParserFast((sample_name, self.mg_no_pmoA_hits_fasta_path),
                                     self.__class__.ref_index,
                                     self.blast_out_no_hits_path, TEST_DATA_DIR, outdir_name="")
                self._parse_empty_file(bhp, sample_name)

    def _parse_empty_file(self, bhp, sample_name):
        out_path = bhp.parse_blast_hits()
        parsed_hits_df = pd.read_csv(out_path,
                                     sep='\t').drop("Unnamed: 0", axis=1).sort_values(by=['quid', 'suid'],
                                                                                      ignore_index=True)
        expected_df = self.blast_out_no_hits_expected.drop(['Contig_nt', 'Contig_GC'], axis=1)
        assert_frame_equal(expected_df, parsed_hits_df)
