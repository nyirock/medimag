from unittest import TestCase

from tools.make_sample_file import extract_sample_names_from_filenames


class TestMakeSampleFile(TestCase):
    def test_extract_sample_names_from_filenames_oneSample_paired1(self):
        l1 = ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R_unpaired_bbduk2.fastq.00.0_0.cor.fasta']
        l1_exp = {'m7_1y2': ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R_unpaired_bbduk2.fastq.00.0_0.cor.fasta']}
        l1_actual = extract_sample_names_from_filenames(l1)
        self.assertEqual(l1_actual, l1_exp)

    def test_extract_sample_names_from_filenames_oneSample_paired2(self):
        l1 = ['Li23127_S1_L004_R2_001.fastq.gz',
              'Li23127_S1_L004_R1_001.fastq.gz',
              'Li23127_S1_L003_R2_001.fastq.gz',
              'Li23127_S1_L002_R2_001.fastq.gz',
              'Li23127_S1_L002_R1_001.fastq.gz',
              'Li23127_S1_L001_R2_001.fastq.gz',
              'Li23127_S1_L001_R1_001.fastq.gz']
        l1_exp = {'Li23127_S1'.lower(): sorted(l1)}
        l1_actual = extract_sample_names_from_filenames(l1)
        self.assertEqual(l1_actual, l1_exp)

    def test_extract_sample_names_from_filenames_twoSamples_paired(self):
        l1 = ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta',
              'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta',
              'Li23127_S1_L004_R2_001.fastq.gz',
              'Li23127_S1_L004_R1_001.fastq.gz',]
        l1_exp = l1_exp = {'Li23127_S1'.lower(): ['Li23127_S1_L004_R1_001.fastq.gz', 'Li23127_S1_L004_R2_001.fastq.gz',],
                           'm7_1y2': ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta']}
        l1_actual = extract_sample_names_from_filenames(l1)
        self.assertEqual(l1_actual, l1_exp)

