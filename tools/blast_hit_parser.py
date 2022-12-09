"""
Placeholder for the hit parser class
Could have a parent HitParser, with a child BlastHitParser to implement mg_wrapser functionality
"""
import json
import os

import numpy as np
from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import GC

from config import ALEN_BP, IDEN, E_VAL, ALEN_PERCENT
from tools.fasta_utils import parse_contigs_ind, get_nt_size, retrive_sequence, write_recruited_reads_to_fasta
from tools.fs import write_df_to_csv, write_dct_to_csv


class BlastHitParser():
    #TODO: maye use inheritance to remove duplication in out_dir_name (inherit both this and BlastRunner from a parent)
    """
    Much slower than a fast parser, but provides sequence data for the recruited reads and GC values
    """
    def __init__(self, sname2fasta_path_tpl, ref_index, raw_blast_hit_df_path, workdir_path,
                 outdir_name="blast", logs_dir_name = "logs",
                 custom_blast_columns_lst=[], genome_stats_log=False):
        self.blast_columns = ['quid', 'suid', 'iden', 'alen',
                              'mism', 'gapo', 'qsta', 'qend',
                              'ssta', 'send', 'eval', 'bits'] + custom_blast_columns_lst

        self.fasta_stats_log = genome_stats_log
        self.logs_dir_name = logs_dir_name
        self.workdir_path = workdir_path
        self.sname, self.input_fasta_path = sname2fasta_path_tpl
        self.ref_index = ref_index
        self.raw_blast_hit_df = self.load_blast_hits_to_df(raw_blast_hit_df_path)
        self.outdir_name = outdir_name
        # print(len(self.ref_index))
        # print(len(self.raw_blast_hit_df))

    def get_sample_name(self):
        return self.sname

    def get_out_dir_path(self):
        return os.path.join(self.workdir_path, self.outdir_name)

    def parse_blast_hits(self):
        input_file_index = parse_contigs_ind(self.input_fasta_path)
        input_file_nt_size = get_nt_size(input_file_index)
        input_file_read_cnt = len(input_file_index)
        #TODO: this is a mistake, it's not a function but a variable: Move log writing to blast_runner, maybe in a separate thread
        if self.write_genome_stats_log:
            self.write_genome_stats_log(input_file_read_cnt, input_file_nt_size)
        recruited_mg = self._unique_queries_top_scored(self.raw_blast_hit_df)
        #recruited_mg = self._unique_queries_top_scored_pandas(self.raw_blast_hit_df)

        self._calculate_metrics(input_file_index, input_file_nt_size, input_file_read_cnt, recruited_mg)

        recruited_mg = self.order_columns(recruited_mg)
        recruited_mg = self.filter_hits(recruited_mg)



        outfile_path = os.path.join(self.workdir_path, self.outdir_name, self.sname)
        outfile_path_csv = outfile_path + ".csv"
        write_df_to_csv(recruited_mg, outfile_path_csv)
        write_recruited_reads_to_fasta(recruited_mg, input_file_index, outfile_path+".fasta")
        input_file_index.close()
        
        return outfile_path_csv

    def write_genome_stats_log(self, n_reads, total_bp, ext=".genome_stats.tsv"):
        dir_path = os.path.join(self.workdir_path, self.outdir_name, self.logs_dir_name)

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        data = {"number_of_reads": n_reads,
                "total_base_pairs": total_bp}
        write_dct_to_csv(data, dir_path, self.sname+ext)


    def load_blast_hits_to_df(self, filepath):
        #these come accordingly with blast runner

        blast_cols = self.blast_columns
        try:
            df = pd.read_csv(filepath, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=blast_cols)  # making an empty dataframe
        df.columns = blast_cols
        return df


    def filter_hits(self, recruited_mg):
        if ALEN_PERCENT:
            recruited_mg = recruited_mg[
                (recruited_mg['iden'] >= IDEN) & (recruited_mg['Coverage'] >= ALEN_PERCENT) & (
                        recruited_mg['eval'] <= E_VAL)]
        elif ALEN_BP:
            recruited_mg = recruited_mg[
                (recruited_mg['iden'] >= IDEN) & (recruited_mg['alen'] >= ALEN_BP) & (
                        recruited_mg['eval'] <= E_VAL)]
        return recruited_mg

    def order_columns(self, recruited_mg):
        return recruited_mg[
            ['quid', 'suid', 'iden', 'alen', 'Coverage', 'Metric', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send',
             'eval', 'bits', 'Ref_size', 'Ref_GC', 'Contig_size', 'Contig_GC', 'Contig_percent', 'Read_RPKM', 'Read_RPKM_wt',
             'Contig_nt']]

    def _calculate_metrics(self, input_file_index, input_file_nt_size, input_file_read_cnt, recruited_mg):
        contig_list = recruited_mg['quid'].tolist()
        recruited_mg['Contig_nt'] = retrive_sequence(contig_list, input_file_index) #nucleotide sequence for the retrieved contig
        recruited_mg['Contig_size'] = recruited_mg['Contig_nt'].apply(lambda x: len(x))
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
        try:
            recruited_mg['Contig_GC'] = recruited_mg['Contig_nt'].apply(lambda x: GC(x))
        except:
            recruited_mg['Contig_GC'] = recruited_mg['Contig_nt'].apply(lambda x: None)

        self._insert_contig_percent_col(input_file_read_cnt, recruited_mg)

        try:
            recruited_mg['Read_RPKM'] = 1.0 / (
                    (recruited_mg['Ref_size'] / 1000.0) * (len(input_file_index) / 1000000.0))
            # weighted Read_RPKM = Read_RPKM * alen
            recruited_mg['Read_RPKM_wt'] = recruited_mg['alen'] / (
                    (recruited_mg['Ref_size'] / 1000.0) * (input_file_nt_size / 1000000.0))
        except:
            recruited_mg['Read_RPKM'] = np.nan

    def _insert_contig_percent_col(self, input_file_read_cnt, recruited_mg):
        assert input_file_read_cnt > 0, "Number of genomic reads has to be more than zero"
        recruited_mg['Contig_percent'] = 100.0 / input_file_read_cnt

    def _unique_queries_top_scored(self, dataframe):
        """returns pandas series object
        each metagenomic read is added once in the order of
         the highest bitscore, lowest eval, highest alen, highest iden.
         In case of ties by all four parameters the first hit is returned.
         ~4x more efficient than using pandas capabilities
         """

        variables = list(dataframe.columns.values)
        scaffolds = dict()
        rows = list()

        for row in dataframe.itertuples():

            # if row[1]=='Ga0073928_10002560':
            if row[1] not in scaffolds:
                #first entry
                scaffolds[row[1]] = row
            else:
                #entering the sorting stage
                #max bits
                if row[12] > scaffolds[row[1]][12]: #12th element is bitscore, as it starts from index
                    scaffolds[row[1]] = row
                #equal bits but min eval
                elif (row[12] == scaffolds[row[1]][12]) and (row[11] < scaffolds[row[1]][11]):
                    scaffolds[row[1]] = row
                elif (row[12] == scaffolds[row[1]][12]) and (row[11] == scaffolds[row[1]][11]) and (row[4] > scaffolds[row[1]][4]):
                    scaffolds[row[1]] = row
                elif (row[12] == scaffolds[row[1]][12]) and (row[11] == scaffolds[row[1]][11]) and (row[4] > scaffolds[row[1]][4]) and (row[3] > scaffolds[row[1]][3]):
                    scaffolds[row[1]] = row
        rows = scaffolds.values()
        df = pd.DataFrame([[getattr(i, j) for j in variables] for i in rows], columns=self.blast_columns)
        return df

    def _unique_queries_top_scored_pandas(self, dataframe):
    # def _unique_queries_top_scored(self, dataframe):
        """Another helper function to decide the best hit within a group
        Sorting a group by ['bits', 'eval', 'alen', 'iden'] gives so far the max amt of hits
        with the default settings
        ~4x slower than the iterative version"""
        variables = list(dataframe.columns.values)
        scaffolds = dict()
        rows = list()


        for i, group in dataframe.groupby('quid'):
            #srt_group = group.sort_values('bits', ascending=False)
            # for row in srt_group.itertuples():
            #     scaffolds[row[1]] = row
            #     break #get the first element and get out

            rows.append(self.__get_best_hit(group))

        #rows = scaffolds.values()
        df = pd.DataFrame(rows)
        return df

    def __get_best_hit(self, df):
        """
        In case of a tie, the first reference sequence is returned
        Gives more sequences than default
        :param df: hits dataframe for a given metagenomic read
        :return: The top hit of input dataframe
        """
        assert len(df) > 0, "Zero length hits group dataframe provided"
        if len(df) == 1:
            return df.iloc[0, :]
        else:
            return df.sort_values(['bits', 'eval', 'alen', 'iden'],
                                      ascending=[False, True, False, False]).iloc[0, :]

