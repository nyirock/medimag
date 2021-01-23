"""
Placeholder for the hit parser class
Could have a parent HitParser, with a child BlastHitParser to implement mg_wrapser functionality
"""
import os

import numpy as np
from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import GC

from config import ALEN_BP, IDEN, E_VAL, ALEN_PERCENT
from tools.fasta_utils import parse_contigs_ind, get_nt_size, retrive_sequence, write_recruited_reads_to_fasta
from tools.fs import write_df_to_csv


class BlastHitParser():
    def __init__(self, sname2fname_tpl, ref_index, raw_blast_hit_df, outdir):
        self.sname, self.input_file_fasta = sname2fname_tpl
        self.ref_index = ref_index
        self.raw_blast_hit_df = raw_blast_hit_df
        self.outdir = outdir
        # print(len(self.ref_index))
        # print(len(self.raw_blast_hit_df))

    def parse_blast_hits(self):
        input_file_index = parse_contigs_ind(self.input_file_fasta)
        input_file_nt_size = get_nt_size(input_file_index)

        recruited_mg = self._unique_scaffold_topBits(self.raw_blast_hit_df)

        self._calculate_metrics(input_file_index, input_file_nt_size, recruited_mg)

        recruited_mg = self.order_columns(recruited_mg)
        recruited_mg = self.filter_hits(recruited_mg)



        outfile = os.path.join(self.outdir, self.sname)
        write_df_to_csv(recruited_mg, outfile+".csv")
        write_recruited_reads_to_fasta(recruited_mg, input_file_index, outfile+".fasta")
        input_file_index.close()


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
             'eval', 'bits', 'Ref_size', 'Ref_GC', 'Contig_size', 'Contig_GC', 'Read_RPKM', 'Read_RPKM_wt',
             'Contig_nt']]

    def _calculate_metrics(self, input_file_index, input_file_nt_size, recruited_mg):
        contig_list = recruited_mg['quid'].tolist()
        recruited_mg['Contig_nt'] = retrive_sequence(contig_list, input_file_index)
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
        try:
            recruited_mg['Read_RPKM'] = 1.0 / (
                    (recruited_mg['Ref_size'] / 1000.0) * (len(input_file_index) / 1000000.0))
            # weighted Read_RPKM = Read_RPKM * alen
            recruited_mg['Read_RPKM_wt'] = recruited_mg['alen'] / (
                    (recruited_mg['Ref_size'] / 1000.0) * (input_file_nt_size / 1000000.0))
        except:
            recruited_mg['Read_RPKM'] = np.nan

    def _unique_scaffold_topBits(self, dataframe):
        # returns pandas series object
        variables = list(dataframe.columns.values)
        scaffolds = dict()
        rows = list()

        for row in dataframe.itertuples():

            # if row[1]=='Ga0073928_10002560':
            if row[1] not in scaffolds:
                scaffolds[row[1]] = row
            else:
                if row[12] > scaffolds[row[1]][12]:
                    scaffolds[row[1]] = row
        rows = scaffolds.values()
        # variables=['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'e_val', 'bits']
        df = pd.DataFrame([[getattr(i, j) for j in variables] for i in rows], columns=variables)
        return df






