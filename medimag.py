"""
The main control over the functional pieces.
Governs blast_runner including
-
"""
import os

#from tools.blast_hit_parser import BlastHitParser

import time
from pprint import pprint

from tools.blast_hit_parser import BlastHitParser
from tools.blast_hit_parser_fast import BlastHitParserFast
from tools.blast_runner import BlastRunner
from config import BLAST_DIR_NAME, BLAST_DB_NAME, WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH
from tools.blast_runner_custom import BlastRunnerCustom
from tools.fastq_reformatter import FastqReformatter
from tools.fs import init_dir, write_sample_file_to_disk, load_sample_file_from_disk, load_blast_hits_to_df, \
    write_dct_to_csv, write_sample_file_dct_to_disk, load_sample_from_disk_to_dict
from tools.fasta_utils import parse_contigs_ind
from tools.make_sample_file import make_sample_file, get_sample_file_dct

# WORKDIR_PATH = os.path.abspath("out")
# INPUT_DIR_PATH = os.path.abspath("/media/andriy/linuxData/mg/m7_mg/_error_corrected_reads/")
WORKDIR_PATH = os.path.abspath("out_m7_basic_workflow2")
#WORKDIR_PATH = os.path.abspath("/media/andriy/SeagateExpansionDrive/out_large_2pp_closest_delta_pmoA/")
INPUT_DIR_PATH = os.path.abspath("in")

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
#WORKDIR_PATH    = os.path.abspath("out_medimag")
# WORKDIR_PATH = "out_medimag_all_m7"
# INPUT_DIR_PATH = os.path.abspath("/media/andriy/linuxData/mg/m7_mg/_error_corrected_reads/")
#########Sample Files###################################################
# sample_file = make_sample_file(INPUT_DIR_PATH, "*fastq.gz")
#write_sample_file_to_disk(sample_file, WORKDIR_PATH)
# sample_file_upd = load_sample_file_from_disk(WORKDIR_PATH)
# assert sample_file.equals(sample_file_upd), "testing if loading is done correctly"
#sample_file_dct = get_sample_file_dct(INPUT_DIR_PATH, "*fastq*.gz")

def fastq_reformat_workflow(sample_file_dct,
                            workdir_path,
                            input_dir_path):
    fr = FastqReformatter(workdir_path)

    # file_path_lst = list(map(lambda x: os.path.join(INPUT_DIR_PATH, x), file_lst))
    # print(file_path_lst)
    # to_process = {"m7_1y":  ["m7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz",
    #                          "m7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz",
    #                          "m7_1y_R_unpaired.00.0_0.cor.fastq.gz"]}
    # print(to_process)

    processed_fastq = {}
    for sname, file_lst in sample_file_dct.items():
        processed_fastq[sname] = fr.run_fastq_reformatter(
            list(map(lambda x: os.path.join(input_dir_path, x), file_lst)), sname)

    pprint(processed_fastq)
    write_dct_to_csv(processed_fastq, fr.out_fasta_dir_path, ".fastq_reformat.log")
    return processed_fastq

def custom_blast_workflow(sname_to_fasta_dict, reference_path, workdir_path):
    blast_runner = BlastRunnerCustom(reference_path, workdir_path, custom_headers="qlen")
    blasted_dct = {}
    for sname, file_path in sname_to_fasta_dict.items():
        blasted_dct[sname] = blast_runner.run_blast_parallel(file_path, sname)

    pprint(blasted_dct)
    write_dct_to_csv(blasted_dct, blast_runner.blast_dir_path, ".blast_runner.log")
    return blasted_dct

def blast_workflow(sname_to_fasta_dict, reference_path, workdir_path):
    blast_runner = BlastRunner(reference_path, workdir_path)
    blasted_dct = {}
    for sname, file_path in sname_to_fasta_dict.items():
        blasted_dct[sname] = blast_runner.run_blast_parallel(file_path, sname)

    pprint(blasted_dct)
    write_dct_to_csv(blasted_dct, blast_runner.blast_dir_path, ".blast_runner.log")
    return blasted_dct

def parse_blast_hits_workflow(sname_to_fasta_dict, blasted_dct, reference_path, workdir_path):
    print("Parsing Blast Hits Started")
    start = time.time()
    ref_index = parse_contigs_ind(reference_path)
    hit_parsers = []
    for sname, fasta_path in sname_to_fasta_dict.items():
        hit_parsers.append(BlastHitParser((sname, fasta_path),
                                          ref_index,
                                          blasted_dct[sname],
                                          workdir_path))
    parsed_dct = {}
    for hp in hit_parsers:
        parsed_dct[hp.get_sample_name()] = hp.parse_blast_hits()
    print(f"Hits Parsing took: {time.time() - start}s")

    pprint(parsed_dct)
    write_dct_to_csv(parsed_dct, hp.get_out_dir_path(), ".hit_parser.log")
    return parsed_dct

def parse_blast_hits_fast_workflow(sname_to_fasta_dict, blasted_dct, reference_path, workdir_path):
    print("Parsing Blast Hits Started")
    start = time.time()
    ref_index = parse_contigs_ind(reference_path)
    hit_parsers = []
    for sname, fasta_path in sname_to_fasta_dict.items():
        hit_parsers.append(BlastHitParserFast((sname, fasta_path),
                                              ref_index,
                                              blasted_dct[sname],
                                              workdir_path,
                                              custom_blast_columns_lst=['qlen']))
    parsed_dct = {}
    for hp in hit_parsers:
        parsed_dct[hp.get_sample_name()] = hp.parse_blast_hits()
    print(f"Hits Parsing took: {time.time() - start}s")

    pprint(parsed_dct)
    write_dct_to_csv(parsed_dct, hp.get_out_dir_path(), ".hit_parser.log")
    return parsed_dct

def init_sample_file(workdir_path, input_dir_path, file_mask="*", identify_pairs=True,
                     sample_file_name="sample_file.tsv", write_to_disk=True):
    init_dir(workdir_path)
    sample_file_dct = get_sample_file_dct(input_dir_path, ext=file_mask, identify_pairs=identify_pairs)
    if write_to_disk:
        write_sample_file_dct_to_disk(sample_file_dct, workdir_path, file_name=sample_file_name)
    return sample_file_dct

##################################################################################
####################High-Level Workflows##########################################
##################################################################################

def run_basic_workflow(workdir_path, input_dir_path, reference_path, sample_file_dict=None,
                       identify_pairs=False, file_mask="*fastq*gz"):
    if not sample_file_dict:
        sf = init_sample_file(workdir_path, input_dir_path, file_mask=file_mask, identify_pairs=identify_pairs)
    else:
        sf = sample_file_dict
    pprint(sf)
    refomat_dct = fastq_reformat_workflow(sf, workdir_path, input_dir_path)
    # reformat_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, ".rf_samples"),
    #                                "*rf.fasta", identify_pairs=False,
    #                               write_to_disk=False)
    # reformat_sf2 = add_path_to_sample_dct(reformat_sf, os.path.join(WORKDIR_PATH, ".rf_samples"))
    # pprint(reformat_sf2)
    # pprint(refomat_dct)
    blasted_dct = blast_workflow(refomat_dct, reference_path, workdir_path)
    pprint(blasted_dct)
    # blasted_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, "blast"), "*tab",
    #                               identify_pairs=False, write_to_disk=False)
    parsed_hits_dct = parse_blast_hits_workflow(refomat_dct, blasted_dct, reference_path, workdir_path)
    pprint(parsed_hits_dct)


def run_fast_workflow(workdir_path, input_dir_path, reference_path,
                      sample_file_dict=None, identify_pairs=False, file_mask="*fastq*gz"):
    #init sample file from dictionary if none was supplied
    if not sample_file_dict:
        sf = init_sample_file(workdir_path, input_dir_path, file_mask=file_mask, identify_pairs=identify_pairs)
    else:
        sf = sample_file_dict
    pprint(sf)
    refomat_dct = fastq_reformat_workflow(sf, workdir_path, input_dir_path)
    # reformat_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, ".rf_samples"),
    #                                "*rf.fasta", identify_pairs=False,
    #                               write_to_disk=False)
    # reformat_sf2 = add_path_to_sample_dct(reformat_sf, os.path.join(WORKDIR_PATH, ".rf_samples"))
    # pprint(reformat_sf2)
    # pprint(refomat_dct)
    blasted_dct = custom_blast_workflow(refomat_dct, reference_path, workdir_path)
    pprint(blasted_dct)
    # blasted_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, "blast"), "*tab",
    #                               identify_pairs=False, write_to_disk=False)
    parsed_hits_dct = parse_blast_hits_fast_workflow(refomat_dct, blasted_dct, reference_path, workdir_path)
    pprint(parsed_hits_dct)

# def main():
#     ########################################################################
#     sample_file_dct = load_sample_from_disk_to_dict(WORKDIR_PATH)
#     pprint(sample_file_dct)
#
#     #########Running Fasta Reformat#########################################
#
#     fr = FastqReformatter(WORKDIR_PATH)
#
#     #file_path_lst = list(map(lambda x: os.path.join(INPUT_DIR_PATH, x), file_lst))
#     #print(file_path_lst)
#     # to_process = {"m7_1y":  ["m7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz",
#     #                          "m7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz",
#     #                          "m7_1y_R_unpaired.00.0_0.cor.fastq.gz"]}
#     # print(to_process)
#
#     processed_fastq = {}
#     for sname, file_lst in sample_file_dct.items():
#          processed_fastq[sname] = fr.run_fastq_reformatter(list(map(lambda x: os.path.join(INPUT_DIR_PATH, x), file_lst)), sname)
#
#     pprint(processed_fastq)
#     write_dct_to_csv(processed_fastq, fr.out_fasta_dir_path, ".fastq_reformat.log")
#     ########################################################################
#
#
#
#
#     #########Running Blast###################################################
#     # sname = "m7_1y2"
#     blast_runner = BlastRunner(REFERENCE_PATH, WORKDIR_PATH)
#     blasted_dct = {}
#     for sname, file_path in processed_fastq.items():
#          blasted_dct[sname] = blast_runner.run_blast_parallel(file_path, sname)
#
#     pprint(blasted_dct)
#     write_dct_to_csv(blasted_dct, blast_runner.blast_dir_path, ".blast_runner.log")
#     ##########################################################################
#
#
#     ##########Running Blast Hit Parser########################################
#     print("Parsing Blast Hits Started")
#     start = time.time()
#     ref_index = parse_contigs_ind(REFERENCE_PATH)
#     hit_parsers = []
#     for sname, fasta_path in processed_fastq.items():
#          hit_parsers.append(BlastHitParser((sname, fasta_path),
#                                            ref_index,
#                                            load_blast_hits_to_df(blasted_dct[sname]),
#                                            WORKDIR_PATH))
#     parsed_dct = {}
#     for hp in hit_parsers:
#         parsed_dct[hp.get_sample_name()] = hp.parse_blast_hits()
#     print(f"Hits Parsing took: {time.time() - start}s")
#
#     pprint(parsed_dct)
#     write_dct_to_csv(parsed_dct, hp.get_out_dir_path(), ".hit_parser.log")
#     # raw_blast_hits_df = load_blast_hits_to_df(raw_blast_hits_path)
#     # sname2processed_fasta_tpl = list(sample_name_to_processed_files_dct.items())[0]
#     # bhp.parse_blast_hits()
#
#     ##########################################################################




    #
    #
    # sample_file_dct = get_sample_file_dct(INPUT_DIR_PATH, "*fast*")
    #
    # #TODO: this is for the testing purposes only
    # sample_name_to_processed_files_dct = {key: os.path.join(INPUT_DIR_PATH, val[0])for key,val in sample_file_dct.items()}
    #
    # print(sample_file_dct)

    ### Parsing blast hits


    # ref_index = parse_contigs_ind(REFERENCE_PATH)
    # raw_blast_hits_df = load_blast_hits_to_df(raw_blast_hits_path)
    # sname2processed_fasta_tpl = list(sample_name_to_processed_files_dct.items())[0]
    # bhp = BlastHitParser(sname2processed_fasta_tpl, ref_index, raw_blast_hits_df, blast_dir)
    # bhp.parse_blast_hits()
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


def add_path_to_sample_dct(sample_dct, dir_path):
    return {key: os.path.join(dir_path, val[0]) for key, val in sample_dct.items()}


    # #main()
    # sf = init_sample_file(WORKDIR_PATH, INPUT_DIR_PATH, "*fastq*gz", identify_pairs=False)
    # pprint(sf)
    # refomat_dct = fastq_reformat_workflow(sf, WORKDIR_PATH, INPUT_DIR_PATH)
    # # reformat_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, ".rf_samples"),
    # #                                "*rf.fasta", identify_pairs=False,
    # #                               write_to_disk=False)
    # # reformat_sf2 = add_path_to_sample_dct(reformat_sf, os.path.join(WORKDIR_PATH, ".rf_samples"))
    # # pprint(reformat_sf2)
    # #pprint(refomat_dct)
    # blasted_dct = blast_workflow(refomat_dct, REFERENCE_PATH, WORKDIR_PATH)
    # pprint(blasted_dct)
    # # blasted_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, "blast"), "*tab",
    # #                               identify_pairs=False, write_to_disk=False)
    # parsed_hits_dct = parse_blast_hits_workflow(refomat_dct, blasted_dct, REFERENCE_PATH, WORKDIR_PATH)
    # pprint(parsed_hits_dct)
    # # blasted_sf2 = add_path_to_sample_dct(blasted_sf, os.path.join(WORKDIR_PATH, "blast"))
    # #pprint(blasted_sf)
    #
    # # pprint(blasted_sf2)
    # #parsed_dct = parse_blast_hits_workflow(reformat_sf2, blasted_sf2, REFERENCE_PATH, WORKDIR_PATH)
    #
    # # init_sample_file("/media/andriy/linuxData/mg/CAB_diversity/medimag/out_medimag/",
    # #                  "/media/andriy/linuxData/mg/CAB_diversity/medimag/out_medimag/blast/",
    # #                  "*csv", identify_pairs=False, sample_file_name="parsed_hits.tsv")


if __name__ == "__main__":
    init_sample_file(WORKDIR_PATH, INPUT_DIR_PATH,identify_pairs=False, file_mask="*fastq.gz")
    sample_file = load_sample_from_disk_to_dict(WORKDIR_PATH)
    pprint(sample_file)
    run_basic_workflow(WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH, sample_file_dict=sample_file)
    #run_fast_workflow(WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH, sample_file_dict=sample_file)