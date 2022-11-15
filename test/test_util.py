from unittest import TestCase

from tools.utils import remove_common_substring

cases = [("123.fq.gz", "123"), ("1234.fq.fastq.gz", "1234"), ("1234.fastq.fq.gz", "1234")]

class TestUtil(TestCase):
    def test_remove_common_suffix(self):
        for filename, expected in cases:
            with self.subTest(msg=filename):
                self.assertEqual(expected, remove_common_substring(filename))

