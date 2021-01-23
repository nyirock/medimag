"""
The main control over the functional pieces.
Governs blast_runner including
-
"""
import os

#from tools.blast_hit_parser import BlastHitParser
from tools.blast_hit_parser import BlastHitParser
from tools.blast_runner import BlastRunner
from config import BLAST_DIR_NAME, BLAST_DB_NAME, WORKDIR, INPUT_DIR, REFERENCE
from tools.fs import init_blast_dir, write_sample_file_to_disk, load_sample_file_from_disk, load_blast_hits_to_df
from tools.fasta_utils import parse_contigs_ind
from tools.make_sample_file import make_sample_file, get_sample_file_dct



#input_files = glob("in/*.csv")

# sample2filename = {'m71y2': ['in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_0.csv',
#                              'in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_1.csv',
#                              'in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_2.csv',],
#                    'm71y': ['in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_3.csv',
#                              'in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_4.csv',
#                              'in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_5.csv',],
#                    'm7': ['in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_6.csv',
#                              'in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_7.csv',
#                              'in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_8.csv',],
#                    }

#print(input_files)

# m71y2_diversity = DiversityCalculator(sample2filename['m71y2'], 'm71y2')
# df_0_diversity = DiversityCalculator(['in/m71y2_m71y_m7_to_pmoA_e_0.0001_iden_45.0%_alen_70.0%_recruited_mg_0.csv'],
#                                      'df_0')

blast_dir = raw_blast_hits_path = os.path.join(WORKDIR, BLAST_DIR_NAME)
raw_blast_hits_path = os.path.join(WORKDIR, BLAST_DIR_NAME, "m7_1y2.tab")

db_path = init_blast_dir(WORKDIR, BLAST_DIR_NAME, BLAST_DB_NAME)
#db_path = os.path.join(WORKDIR, BLAST_DIR_NAME, BLAST_DB_NAME)
blast_runner = BlastRunner(REFERENCE, db_path)
blast_runner.run_blast_parallel(os.path.join(INPUT_DIR, "m7_1y2_R1_bbduk2.fastq.00.0_0.cor.fasta"),
                                raw_blast_hits_path)
# print(os.path.join(WORKDIR, BLAST_DIR_NAME, BLAST_DB_NAME))

sample_file = make_sample_file(INPUT_DIR, "*fastq.gz")
write_sample_file_to_disk(sample_file, WORKDIR)
sample_file_upd = load_sample_file_from_disk(WORKDIR)
assert sample_file.equals(sample_file_upd), "testing if loading is done correctly"


sample_file_dct = get_sample_file_dct(INPUT_DIR, "*fast*")

#TODO: this is for the testing purposes only
sample_name_to_processed_files_dct = {key: os.path.join(INPUT_DIR, val[0])for key,val in sample_file_dct.items()}


### Parsing blast hits
print(sample_file_dct)
ref_index = parse_contigs_ind(REFERENCE)
raw_blast_hits_df = load_blast_hits_to_df(raw_blast_hits_path)
sname2processed_fasta_tpl = list(sample_name_to_processed_files_dct.items())[0]
bhp = BlastHitParser(sname2processed_fasta_tpl, ref_index, raw_blast_hits_df, blast_dir)
bhp.parse_blast_hits()
# print(df_0_diversity.concordance)
#
# print(m71y2_diversity.concordance)
#
# diversity_lst = []
# for sample_name, files_lst in sample2filename.items():
#     dc = DiversityCalculator(files_lst, sample_name)
#     if dc.no_hits:
#         continue
#     diversity_lst.append(dc)
#
# sorted_lst = sorted(diversity_lst, key=lambda x: x.concordance.mean)
#
# print(sorted_lst)
#
# for i in sorted_lst:
#     print(i)