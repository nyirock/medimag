"""
Fast version of blast hit parser that doesn't require making index of the input file.
REQUIRES: qlen blast header is added to blast headers:
"""
import os

import numpy as np
from Bio.SeqUtils import GC

from tools.blast_hit_parser import BlastHitParser
from tools.fasta_utils import get_nt_size_ext, run_seqkit_stats, get_nt_size_read_cnt_from_seqkit
from tools.fs import write_df_to_csv


class BlastHitParserFast(BlastHitParser):
    def __init__(self, *args, custom_blast_columns_lst=["qlen"], **kwargs):
        #self.blast_columns = self.blast_columns + custom_blast_columns_lst
        super(BlastHitParserFast, self).__init__(*args,
                                                 custom_blast_columns_lst=custom_blast_columns_lst,
                                                 **kwargs)

    #TODO: overload the critical bits of code to avoid making an index
    def parse_blast_hits(self):
        """
        No fasta output for the recruited reads
        :return:
        """
        #input_file_index = parse_contigs_ind(self.input_fasta_path)

        input_file_read_cnt, input_file_nt_size = get_nt_size_ext(self.input_fasta_path)
        if self.write_genome_stats_log:
            self.write_genome_stats_log(input_file_read_cnt, input_file_nt_size)


        recruited_mg = self._unique_queries_top_scored(self.raw_blast_hit_df)

        self._calculate_metrics(input_file_read_cnt, input_file_nt_size, recruited_mg)

        recruited_mg = self.order_columns(recruited_mg)
        recruited_mg = self.filter_hits(recruited_mg)

        outfile_path = os.path.join(self.workdir_path, self.outdir_name, self.sname)
        outfile_path_csv = outfile_path + ".csv"
        write_df_to_csv(recruited_mg, outfile_path_csv) 


        return outfile_path_csv

    def total_base_pairs_calculate_metrics(self, input_file_read_cnt, input_file_nt_size, recruited_mg):
        """
        Two columns 'Contig_nt', 'Contig_GC'  are removed here compared with the parent method
        :param input_file_index:
        :param input_file_nt_size:
        :param recruited_mg:
        :return:
        """
        #TODO: add specific exceptions instead of general ones
        contig_list = recruited_mg['quid'].tolist()
        #  ----recruited_mg['Contig_nt'] = retrive_sequence(contig_list, input_file_index)
        # Renaming 'qlen' to 'Contig_size'
        recruited_mg.rename(columns={'qlen': 'Contig_size'}, inplace=True)
        #recruited_mg['Contig_size'] = recruited_mg['qlen']
        # recruited_mg[i]['Ref_nt']=recruited_mg[i]['suid'].apply(lambda x: self.ref_index[str(x)].seq)
        recruited_mg['Ref_size'] = recruited_mg['suid'].apply(lambda x: len(self.ref_index[str(x)]))
        recruited_mg['Ref_GC'] = recruited_mg['suid'].apply(lambda x: GC(self.ref_index[str(x)].seq))
        # recruited_mg[i]['Coverage']=recruited_mg[i]['alen'].apply(lambda x: 100.0*float(x))/min(recruited_mg[i]['Contig_size'].apply(lambda y: y),recruited_mg[i]['Ref_size'].apply(lambda z: z))
        # df.loc[:, ['B0', 'B1', 'B2']].min(axis=1)
        recruited_mg['Coverage'] = recruited_mg['alen'].apply(lambda x: 100.0 * float(x)) / recruited_mg.loc[:,
                                                                                            ["Contig_size",
                                                                                             "Ref_size"]].min(
            axis=1)
        recruited_mg['Metric'] = recruited_mg['Coverage'] * recruited_mg['iden'] / 100.0
        #----- try:
        #-----     recruited_mg['Contig_GC'] = recruited_mg['Contig_nt'].apply(lambda x: GC(x))
        #----- except:
        #-----     recruited_mg['Contig_GC'] = recruited_mg['Contig_nt'].apply(lambda x: None)
        try:
            recruited_mg['Read_RPKM'] = 1.0 / (
                    (recruited_mg['Ref_size'] / 1000.0) * (input_file_read_cnt / 1000000.0))
            # weighted Read_RPKM = Read_RPKM * alen
            recruited_mg['Read_RPKM_wt'] = recruited_mg['alen'] / (
                    (recruited_mg['Ref_size'] / 1000.0) * (input_file_nt_size / 1000000.0))
        except:
            recruited_mg['Read_RPKM'] = np.nan

    def order_columns(self, recruited_mg):
        return recruited_mg[
            ['quid', 'suid', 'iden', 'alen', 'Coverage', 'Metric', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send',
             'eval', 'bits', 'Ref_size', 'Ref_GC', 'Contig_size', 'Read_RPKM', 'Read_RPKM_wt']]