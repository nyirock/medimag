"""
The main control over the functional pieces.
Governs blast_runner including
-
"""
import os

#from tools.blast_hit_parser import BlastHitParser

import time
from collections import defaultdict
from pprint import pprint
import pandas as pd

from tools.blast_hit_parser import BlastHitParser
from tools.blast_hit_parser_fast import BlastHitParserFast
from tools.blast_runner import BlastRunner
from config import BLAST_DIR_NAME, BLAST_DB_NAME, WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH
from tools.blast_runner_custom import BlastRunnerCustom
from tools.diversity_calculator import DiversityCalculator
from tools.fastq_qc import FastqQC
from tools.fastq_reformatter import FastqReformatter
from tools.fs import init_dir, write_sample_file_to_disk, load_sample_file_from_disk, load_blast_hits_to_df, \
    write_dct_to_csv, write_sample_file_dct_to_disk, load_sample_from_disk_to_dict, write_df_to_csv
from tools.fasta_utils import parse_contigs_ind
from tools.make_sample_file import make_sample_file, get_sample_file_dct

#WORKDIR_PATH = os.path.abspath("/media/andriy/SeagateExpansionDrive/out_medimag_pp_old_and_new_nmd_forest_river_qc")
WORKDIR_PATH = os.path.abspath("out_16s_mg_large_fast")
#INPUT_DIR_PATH = os.path.abspath("/media/andriy/SeagateExpansionDrive/pp_old_and_new_nmd_forest_river_qc/")
INPUT_DIR_PATH = os.path.abspath("in_mg_large/")
# INPUT_DIR_PATH = os.path.abspath("/media/andriy/linuxData/mg/m7_mg/_error_corrected_reads/")
#WORKDIR_PATH = os.path.abspath("several_mg_raw_doubles_combined_basic_w_qc_no_gz_paired")
#WORKDIR_PATH = os.path.abspath("/media/andriy/SeagateExpansionDrive/out_large_2pp_closest_delta_pmoA/")
#INPUT_DIR_PATH = os.path.abspath("several_mg_raw_basic_w_qc_gz_paired/.qc_samples/")

# INPUT_DIR_PATH="in_gz_glitch"
# WORKDIR_PATH="in_gz_glitch"
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

def fastq_qc_workflow(sample_file_dct,
                            workdir_path,
                            input_dir_path,
                            gzipped=True):
    fqc = FastqQC(workdir_path, gzipped=gzipped)

    # file_path_lst = list(map(lambda x: os.path.join(INPUT_DIR_PATH, x), file_lst))
    # print(file_path_lst)
    # to_process = {"m7_1y":  ["m7_1y_R1_bbduk2.fastq.00.0_0.cor.fastq.gz",
    #                          "m7_1y_R2_bbduk2.fastq.00.0_0.cor.fastq.gz",
    #                          "m7_1y_R_unpaired.00.0_0.cor.fastq.gz"]}
    # print(to_process)

    qc_fastq = defaultdict(list)
    for sname, file_lst in sample_file_dct.items():
        assert len(file_lst) <= 2, "Max of 2 files are allowed per sample"
        for sname, file_lst in sample_file_dct.items():
            qc_fastq[sname] = fqc.run_fastq_qc(
                list(map(lambda x: os.path.join(input_dir_path, x), file_lst)), sname)
        pprint(qc_fastq)
        write_dct_to_csv(qc_fastq, fqc.out_fasta_dir_path, ".fastq_qc.log")
        return qc_fastq


def fastq_reformat_workflow(sample_file_dct,
                            workdir_path,
                            input_dir_path,
                            gzipped=False):
    fr = FastqReformatter(workdir_path, gzipped=gzipped)

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

def parse_blast_hits_fast_workflow(sname_to_fasta_dict, blasted_dct, reference_path, workdir_path, genome_stats_log=True):
    print("Parsing Blast Hits Started")
    start = time.time()
    ref_index = parse_contigs_ind(reference_path)
    hit_parsers = []
    for sname, fasta_path in sname_to_fasta_dict.items():
        hit_parsers.append(BlastHitParserFast((sname, fasta_path),
                                              ref_index,
                                              blasted_dct[sname],
                                              workdir_path,
                                              custom_blast_columns_lst=['qlen'],
                                              genome_stats_log= genome_stats_log))
    parsed_dct = {}
    for hp in hit_parsers:
        parsed_dct[hp.get_sample_name()] = hp.parse_blast_hits()
    print(f"Hits Parsing took: {time.time() - start}s")

    pprint(parsed_dct)
    write_dct_to_csv(parsed_dct, hp.get_out_dir_path(), ".hit_parser.log")
    return parsed_dct

def init_sample_file(workdir_path, input_dir_path, file_mask="*", identify_pairs=True, strict_pairs=False,
                     sample_file_name="sample_file.tsv", write_to_disk=True):
    init_dir(workdir_path)
    sample_file_dct = get_sample_file_dct(input_dir_path, ext=file_mask, identify_pairs=identify_pairs, strict_pairs=strict_pairs)
    if write_to_disk:
        write_sample_file_dct_to_disk(sample_file_dct, workdir_path, file_name=sample_file_name)
    return sample_file_dct

##################################################################################
####################High-Level Workflows##########################################
##################################################################################

def run_basic_workflow(workdir_path, input_dir_path, reference_path,
                       sample_file_dict=None, identify_pairs=False, file_mask="*fastq*",
                       gzipped=False, qc=False):
    if not sample_file_dict:
        sf = init_sample_file(workdir_path, input_dir_path, file_mask=file_mask, identify_pairs=identify_pairs)
    else:
        sf = sample_file_dict
    pprint(sf)
    refomat_dct = fastq_reformat_workflow(sf, workdir_path, input_dir_path, gzipped=gzipped)
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


def run_workflow(workdir_path, input_dir_path, reference_path, workflow ="fast",
                 sample_file_dict=None, identify_pairs=False, file_mask="*fastq*", gzipped=False, qc=False, just_qc=False):
    #init sample file from dictionary if none was supplied

    if (not sample_file_dict) and qc:
        sf = init_sample_file(workdir_path, input_dir_path, file_mask=file_mask, identify_pairs=identify_pairs, strict_pairs=True)
    elif (not sample_file_dict) and (not qc):
        sf = init_sample_file(workdir_path, input_dir_path, file_mask=file_mask, identify_pairs=identify_pairs,
                              strict_pairs=False)
    else:
        sf = sample_file_dict
    pprint(sf)
    if qc:
        qc_dct = fastq_qc_workflow(sf, workdir_path, input_dir_path, gzipped=gzipped)
        pprint(qc_dct)
        if just_qc:
            return
        #TODO: use dynamic qc folder instead of a constant "qc_samples"
        #TODO: find a better way for file mask
        if not gzipped:
            file_mask = "*fastq"
        qc_dct2 = init_sample_file(workdir_path, os.path.join(workdir_path, ".qc_samples"), file_mask=file_mask, identify_pairs=identify_pairs, strict_pairs=False)
        pprint(qc_dct2)
        refomat_dct = fastq_reformat_workflow(qc_dct2, workdir_path, os.path.join(workdir_path, ".qc_samples"), gzipped=gzipped)
    else:
        refomat_dct = fastq_reformat_workflow(sf, workdir_path, input_dir_path, gzipped=gzipped)
    # reformat_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, ".rf_samples"),
    #                                "*rf.fasta", identify_pairs=False,
    #                               write_to_disk=False)
    # reformat_sf2 = add_path_to_sample_dct(reformat_sf, os.path.join(WORKDIR_PATH, ".rf_samples"))
    # pprint(reformat_sf2)
    # pprint(refomat_dct)
    if workflow == "fast":
        blasted_dct = custom_blast_workflow(refomat_dct, reference_path, workdir_path)
        pprint(blasted_dct)
        # blasted_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, "blast"), "*tab",
        #                               identify_pairs=False, write_to_disk=False)
        parsed_hits_dct = parse_blast_hits_fast_workflow(refomat_dct, blasted_dct, reference_path, workdir_path)
    elif workflow == "basic":
        blasted_dct = blast_workflow(refomat_dct, reference_path, workdir_path)
        pprint(blasted_dct)
        # blasted_sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, "blast"), "*tab",
        #                               identify_pairs=False, write_to_disk=False)
        parsed_hits_dct = parse_blast_hits_workflow(refomat_dct, blasted_dct, reference_path, workdir_path)
    pprint(parsed_hits_dct)

def add_path_to_sample_dct(sample_dct, dir_path):
    return {key: os.path.join(dir_path, val[0]) for key, val in sample_dct.items()}

if __name__ == "__main__":

    #sf = init_sample_file(WORKDIR_PATH, os.path.join(WORKDIR_PATH, ".qc_samples"), identify_pairs=True, file_mask="*fastq.gz", write_to_disk=True, strict_pairs=False)
    #pprint(sf)
    #sf = load_sample_from_disk_to_dict(WORKDIR_PATH)
    #pprint(sf)
    # #
    #run_workflow(WORKDIR_PATH, os.path.join(WORKDIR_PATH, ".qc_samples"), REFERENCE_PATH, sample_file_dict=sf, gzipped=True, workflow='fast', qc=False)
    #run_workflow(WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH, identify_pairs=False, gzipped=False, workflow='fast', qc=False)
    # run_workflow(WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH, identify_pairs=True, gzipped=False, workflow='fast', qc=False)
    #run_workflow(WORKDIR_PATH, os.path.join(WORKDIR_PATH, ".qc_samples"), REFERENCE_PATH, identify_pairs=True, gzipped=True, workflow='fast',
    #             qc=False)

    #run_workflow(WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH, gzipped=False, qc=True, identify_pairs=True)
    # run_workflow(WORKDIR_PATH, INPUT_DIR_PATH, REFERENCE_PATH, gzipped=True, qc=False,
    #                   identify_pairs=False)

    #Testing Diversity Calculator
    # dc = DiversityCalculator(("m7_1y2", "out_m71y_paired_no_qc_gzipped/blast/m7_1y2.csv"),
    #                          "out_m71y_paired_no_qc_gzipped")
    # print(dc)
    # dc.print_reference_hits_summary()

    #################################################################################3
    #Diversity calculator from csv-files
    # csv_dct = init_sample_file(WORKDIR_PATH, INPUT_DIR_PATH, identify_pairs=False, file_mask="*csv", write_to_disk=True,
    #                  strict_pairs=False)
    #
    # csv_dct_w_path = add_path_to_sample_dct(csv_dct, INPUT_DIR_PATH)
    # pprint(csv_dct_w_path)
    # dc_lst = []
    # for sname,filepath in csv_dct_w_path.items():
    #     dc_tmp = DiversityCalculator((sname, filepath), WORKDIR_PATH)
    #     if dc_tmp.no_hits:
    #         continue
    #     dc_lst.append(dc_tmp)
    #
    # cnt_lst = map(DiversityCalculator.get_counts, dc_lst)
    # names_lst = map(lambda x: x.sname, dc_lst)
    # otu_table = pd.DataFrame(cnt_lst, index=names_lst).fillna(0).T
    # write_df_to_csv(otu_table, os.path.join(WORKDIR_PATH, "otu_table.tsv"))
    ##############################################################################
    # Testing FasgtQC module
    # raw_sample_file_strict_pairs = init_sample_file(WORKDIR_PATH, INPUT_DIR_PATH, sample_file_name="qc_samples.tsv",
    #                                                 file_mask="*fastq*gz", identify_pairs=True, strict_pairs=True)
    # pprint(raw_sample_file_strict_pairs)
    # qc_sample_dct = fastq_qc_workflow(raw_sample_file_strict_pairs, WORKDIR_PATH, INPUT_DIR_PATH)
    # pprint(qc_sample_dct)

    #### Fixing gzip capability

    # blasted_dct = blast_workflow(reformat_sf, REFERENCE_PATH, WORKDIR_PATH)


    # #
    # blasted_sf = init_sample_file(WORKDIR_PATH, WORKDIR_PATH, "*tab", identify_pairs=False, write_to_disk=False)
    # pprint(blasted_sf)
    # blasted_sf = add_path_to_sample_dct(blasted_sf, WORKDIR_PATH)
    # print(blasted_sf)
    # #
    # # pprint(blasted_sf)
    # parsed_dct = parse_blast_hits_workflow(reformat_sf, blasted_dct, REFERENCE_PATH, WORKDIR_PATH)


    ############## CURRENT BLASTING WORKFLOW #######################################
    # reformat_sf = init_sample_file(INPUT_DIR_PATH, INPUT_DIR_PATH, "*.fasta*",
    #                                identify_pairs=False, write_to_disk=False)
    # reformat_sf = add_path_to_sample_dct(reformat_sf, INPUT_DIR_PATH)
    # pprint(reformat_sf)
    # blasted_dct = custom_blast_workflow(reformat_sf, REFERENCE_PATH, WORKDIR_PATH)
    # pprint(blasted_dct)
    # parsed_dct = parse_blast_hits_fast_workflow(reformat_sf, blasted_dct, REFERENCE_PATH, WORKDIR_PATH)
    # pprint(parsed_dct)

    #################################################################################
    #Diversity calculator from csv-files
    # csv_dct = init_sample_file(WORKDIR_PATH, INPUT_DIR_PATH, identify_pairs=False, file_mask="*csv", write_to_disk=True,
    #                  strict_pairs=False)
    #
    # csv_dct_w_path = add_path_to_sample_dct(csv_dct, INPUT_DIR_PATH)
    # pprint(csv_dct_w_path)


    ###### Getting Sample stats #####################
    parsed_dct = {'mg_large': '/home/dunfield/Documents/andriy/medimag/out_16s_mg_large_fast/blast/mg_large.csv'}
    dc_lst = []
    for sname,filepath in parsed_dct.items():
        dc_tmp = DiversityCalculator((sname, filepath), WORKDIR_PATH)
        if dc_tmp.no_hits:
            continue
        dc_lst.append(dc_tmp)
    for dc_obj in dc_lst:
        write_df_to_csv(dc_obj.reference_hit_summary, os.path.join(WORKDIR_PATH, dc_obj.sname+"_hit_summary.tsv"))
    ######## End Sample stats #######################

    #OTU table-specific


    # cnt_lst = map(DiversityCalculator.get_counts, dc_lst)
    # print(dc_lst)
    # names_lst = map(lambda x: x.sname, dc_lst)
    # print(list(names_lst))
    # otu_table = pd.DataFrame(cnt_lst, index=names_lst).fillna(0).T
    # write_df_to_csv(otu_table, os.path.join(WORKDIR_PATH, "otu_table.tsv"))