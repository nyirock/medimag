import os
from unittest import TestCase

from tools.fasta_utils import get_nt_size_ext

TEST_DATA_DIR = os.path.abspath('data')

fasta_paths = [(os.path.join(TEST_DATA_DIR, "mg_large.fasta"), (4839683,1089534308)),
               (os.path.join(TEST_DATA_DIR, "mg_no_pmoA_hits.fasta"), (35311 ,9199439)),
               (os.path.join(TEST_DATA_DIR, "all_pmoA_nr99.fasta"), (864, 675924)),
               (os.path.join(TEST_DATA_DIR, "recruited_mg_8.fasta"), (0, 0))]

class TestFS(TestCase):


    def test_get_nt_size_ext(self):
        for filepath, (exp_nt, exp_sz) in fasta_paths:
            nt, sz = get_nt_size_ext(filepath)
            self.assertEqual(nt, exp_nt)
            self.assertEqual(sz, exp_sz)
