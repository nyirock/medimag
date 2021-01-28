from unittest import TestCase

from tools.make_sample_file import extract_sample_names_from_filenames


class TestMakeSampleFile(TestCase):
    def test_extract_sample_names_from_filenames_oneSample_paired1(self):
        l1 = ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R_unpaired_bbduk2.fastq.00.0_0.cor.fasta']
        l1_exp = {'m7_1y2': ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R_unpaired_bbduk2.fastq.00.0_0.cor.fasta']}
        l1_actual = extract_sample_names_from_filenames(l1, identify_pairs=True)
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
        l1_actual = extract_sample_names_from_filenames(l1, identify_pairs=True)
        self.assertEqual(l1_actual, l1_exp)

    def test_extract_sample_names_from_filenames_oneSample_paired2_strict(self):
        l1 = ['Li23127_S1_L004_R2_001.fastq.gz',
              'Li23127_S1_L004_R1_001.fastq.gz',
              'Li23127_S1_L003_R2_001.fastq.gz',
              'Li23127_S1_L003_R1_001.fastq.gz',
              'Li23127_S1_L002_R2_001.fastq.gz',
              'Li23127_S1_L002_R1_001.fastq.gz',
              'Li23127_S1_L001_R2_001.fastq.gz',
              'Li23127_S1_L001_R1_001.fastq.gz']
        l1_exp = {'Li23127_S1_L001'.lower(): ['Li23127_S1_L001_R1_001.fastq.gz', 'Li23127_S1_L001_R2_001.fastq.gz'],
                  'Li23127_S1_L002'.lower(): ['Li23127_S1_L002_R1_001.fastq.gz', 'Li23127_S1_L002_R2_001.fastq.gz'],
                  'Li23127_S1_L003'.lower(): ['Li23127_S1_L003_R1_001.fastq.gz', 'Li23127_S1_L003_R2_001.fastq.gz'],
                  'Li23127_S1_L004'.lower(): ['Li23127_S1_L004_R1_001.fastq.gz', 'Li23127_S1_L004_R2_001.fastq.gz'],}
        l1_actual = extract_sample_names_from_filenames(l1, identify_pairs=True, strict_pairs=True)
        self.assertEqual(dict(l1_actual), l1_exp)

    def test_extract_sample_names_from_filenames_twoSamples_paired(self):
        l1 = ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta',
              'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta',
              'Li23127_S1_L004_R2_001.fastq.gz',
              'Li23127_S1_L004_R1_001.fastq.gz',]
        l1_exp = l1_exp = {'Li23127_S1'.lower(): ['Li23127_S1_L004_R1_001.fastq.gz', 'Li23127_S1_L004_R2_001.fastq.gz',],
                           'm7_1y2': ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta', 'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta']}
        l1_actual = extract_sample_names_from_filenames(l1, identify_pairs=True)
        self.assertEqual(l1_actual, l1_exp)

    def test_extract_sample_names_from_filenames_complexSamples_paired_strict(self):
        samples = ['m7_R2_bbduk2.fastq.00.0_0.cor.fastq.gz', 'm7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                   'Li23128_S2_L003_R2_001.fastq.gz', 'Li23128_S2_L004_R1_001.fastq.gz',
                   'm7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz', 'm7_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                   'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fastq.gz', 'Li23128_S2_L003_R1_001.fastq.gz',
                   'Li23128_S2_L004_R2_001.fastq.gz', 'Li23128_S2_L002_R1_001.fastq.gz',
                   'Li23128_S2_L002_R2_001.fastq.gz', 'Li23128_S2_L001_R2_001.fastq.gz',
                   'm7_1y_R_unpaired.00.0_0.cor.fastq.gz', 'm7_R_unpaired.00.0_0.cor.fastq.gz',
                   'Li23128_S2_L001_R1_001.fastq.gz', 'm7_1y2_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                   'm7_1y2_R_unpaired.00.0_0.cor.fastq.gz']
        samples_dct_expected = {'li23128_s2_l001': ['Li23128_S2_L001_R1_001.fastq.gz',
                                                    'Li23128_S2_L001_R2_001.fastq.gz'],
                                'li23128_s2_l002': ['Li23128_S2_L002_R1_001.fastq.gz',
                                                    'Li23128_S2_L002_R2_001.fastq.gz'],
                                'li23128_s2_l003': ['Li23128_S2_L003_R1_001.fastq.gz',
                                                    'Li23128_S2_L003_R2_001.fastq.gz'],
                                'li23128_s2_l004': ['Li23128_S2_L004_R1_001.fastq.gz',
                                                    'Li23128_S2_L004_R2_001.fastq.gz'],
                                'm7': ['m7_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                       'm7_R2_bbduk2.fastq.00.0_0.cor.fastq.gz'],
                                'm7_1y': ['m7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                          'm7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz'],
                                'm7_1y2': ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                           'm7_1y2_R2_bbduk2.fastq.00.0_0.cor.fastq.gz'],
                                'm7_1y2_r_unpaired_00_0_0_cor': ['m7_1y2_R_unpaired.00.0_0.cor.fastq.gz'],
                                'm7_1y_r_unpaired_00_0_0_cor': ['m7_1y_R_unpaired.00.0_0.cor.fastq.gz'],
                                'm7_r_unpaired_00_0_0_cor': ['m7_R_unpaired.00.0_0.cor.fastq.gz']}

        samples_dct_actual = extract_sample_names_from_filenames(samples, identify_pairs=True, strict_pairs=True)
        self.assertEqual(samples_dct_actual, samples_dct_expected)

    def test_extract_sample_names_from_filenames_complexSamples_paired(self):
        samples = ['m7_R2_bbduk2.fastq.00.0_0.cor.fastq.gz', 'm7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                   'Li23128_S2_L003_R2_001.fastq.gz', 'Li23128_S2_L004_R1_001.fastq.gz',
                   'm7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz', 'm7_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                   'm7_1y2_R1_bbduk2.fastq.00.0_0.cor.fastq.gz', 'Li23128_S2_L003_R1_001.fastq.gz',
                   'Li23128_S2_L004_R2_001.fastq.gz', 'Li23128_S2_L002_R1_001.fastq.gz',
                   'Li23128_S2_L002_R2_001.fastq.gz', 'Li23128_S2_L001_R2_001.fastq.gz',
                   'm7_1y_R_unpaired.00.0_0.cor.fastq.gz', 'm7_R_unpaired.00.0_0.cor.fastq.gz',
                   'Li23128_S2_L001_R1_001.fastq.gz', 'm7_1y2_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                   'm7_1y2_R_unpaired.00.0_0.cor.fastq.gz']
        samples_dct_expected = {'li23128_s2': ['Li23128_S2_L001_R1_001.fastq.gz',
                                               'Li23128_S2_L001_R2_001.fastq.gz',
                                               'Li23128_S2_L002_R1_001.fastq.gz',
                                               'Li23128_S2_L002_R2_001.fastq.gz',
                                               'Li23128_S2_L003_R1_001.fastq.gz',
                                               'Li23128_S2_L003_R2_001.fastq.gz',
                                               'Li23128_S2_L004_R1_001.fastq.gz',
                                               'Li23128_S2_L004_R2_001.fastq.gz'],
                                'm7': ['m7_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                       'm7_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                       'm7_R_unpaired.00.0_0.cor.fastq.gz'],
                                'm7_1y': ['m7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                          'm7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                          'm7_1y_R_unpaired.00.0_0.cor.fastq.gz'],
                                'm7_1y2': ['m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                           'm7_1y2_R2_bbduk2.fastq.00.0_0.cor.fastq.gz',
                                           'm7_1y2_R_unpaired.00.0_0.cor.fastq.gz']}

        samples_dct_actual = extract_sample_names_from_filenames(samples, identify_pairs=True, strict_pairs=False)
        self.assertEqual(samples_dct_actual, samples_dct_expected)



