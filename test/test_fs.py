import os
from unittest import TestCase
import pandas as pd
from pandas._testing import assert_frame_equal

from tools.fs import write_sample_file_to_disk, load_sample_file_from_disk, load_sample_from_disk_to_dict, \
    load_blast_hits_to_df
from tools.make_sample_file import snames2filenames_dict_to_df

TEST_DATA_DIR = 'data'


class TestFS(TestCase):
    def setUp(self) -> None:
        self.snames2files = {'li23127_s1': ['Li23127_S1_L001_R1_001.fastq.gz',
                                            'Li23127_S1_L001_R2_001.fastq.gz',
                                            'Li23127_S1_L002_R1_001.fastq.gz',
                                            'Li23127_S1_L002_R2_001.fastq.gz',
                                            'Li23127_S1_L003_R2_001.fastq.gz',
                                            'Li23127_S1_L004_R1_001.fastq.gz',
                                            'Li23127_S1_L004_R2_001.fastq.gz'],
                             'm7_1y2': ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                        'm7_1y2_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                        'm7_1y2_R_unpaired.00.0_0.cor.fastq.gz'],
                             'm7_1y': ['m7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                       'm7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                       'm7_1y_R_unpaired.00.0_0.cor.fastq.gz'],
                             'm7': ['m7_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                    'm7_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                    'm7_R_unpaired.00.0_0.cor.fastq.gz']}
        self.snames2files_df = snames2filenames_dict_to_df(self.snames2files)

    def test_load_sample_from_disk_to_dict(self):
        write_sample_file_to_disk(self.snames2files_df, TEST_DATA_DIR)
        dct = load_sample_from_disk_to_dict(TEST_DATA_DIR)
        self.assertEqual(self.snames2files, dct)

    def test_load_sample_from_disk(self):
        write_sample_file_to_disk(self.snames2files_df, TEST_DATA_DIR)
        df2 = load_sample_file_from_disk(TEST_DATA_DIR)
        assert_frame_equal(self.snames2files_df, df2)

    def test_load_blast_hits_to_df(self):
        empty_blast_out = os.path.join(TEST_DATA_DIR, "recruited_mg_8.tab")
        empty_blast_out_expected = pd.read_csv(os.path.join(TEST_DATA_DIR, "recruited_mg_8.csv"), sep='\t')
        with self.assertRaises(pd.errors.EmptyDataError):
            pd.read_csv(empty_blast_out)
        df = load_blast_hits_to_df(empty_blast_out)
        assert_frame_equal(pd.DataFrame(columns=['quid', 'suid', 'iden',
                                                 'alen', 'mism', 'gapo',
                                                 'qsta', 'qend', 'ssta',
                                                 'send', 'eval', 'bits']), df)