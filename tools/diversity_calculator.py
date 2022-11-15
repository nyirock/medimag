import os
from collections import namedtuple
from typing import List
import pandas as pd
import numpy as np
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW

from tools.fs import init_dir

Concordance = namedtuple('Concordance', 'mean std sem')

class DiversityCalculator(object):
    """
    Performs parsing of the recruited hits stats, calculating copy numbers in a sample
    with different metrics:
    - Total reads,
    - Total RPKM,
    - Theoretical Copies and,
    - Total Weighted RPKM
    per a database sequence. Saves the results in a feature table.

    Computes diversity metrics for the samples:
    - Concordance: weighted averages on sequence identities by their weighted RPKM
    - Novelty: 1 - concordance (However I think it's better to have it
    as a fraction of the reads that are below average concordance)
    - Target indices: TI, MOXI, AOXI: determined by the phylogenetic distance
    from a target group and then taking a fraction.
    Common diversity indices from scikit-bio like:
    - Shannon, Chao1, Faith_pd, etc. Need to convert to integer counts if needed.
    """
    def __init__(self, sname2parsed_hits_path_tpl, workdir_path, fasta_path=None, outdir_name="diversity"):
        """
        Import csv-files as a dataframe and merge the resultant dataframe.
        :param sname2fasta_path_tpl:  Tuple containing sample name and fasta path: as in fastq_reformat_workflow output
        """
        self.no_hits = False
        self.workdir_path = workdir_path
        self.sname, self.reference_hits_path = sname2parsed_hits_path_tpl
        self.fasta_path = fasta_path
        self.outdir_path = os.path.join(workdir_path, outdir_name)
        self.reference_hits = self._load_hits(self.reference_hits_path)
        if not len(self.reference_hits):
            self.no_hits = True
            #TODO: implement better the case for no hits
            return
        self._convert_iden_percent_to_frac()
        self.reference_hit_summary = self._summarize_reference_hits(self.reference_hits)
        self.concordance = Concordance(*self._calculate_concordance_stats())
        self.novelty = self.calculate_novelty()


    ######################################
    ########### Getters ##################
    ######################################

    def get_counts(self, by='Total_RPKM_wt'):
        assert by in self.reference_hit_summary.columns, "counts_type not recognized"
        return self.reference_hit_summary[by]

    ######################################
    ########### Overloaded methods #######
    ######################################

    def __str__(self):
        return f"Sample name: {self.sname}\n" \
               f"\tConcordance Statistics:\n" \
               f"{'Concordance mean:':<35}{self.concordance.mean:>10}\n" \
               f"{'Concordance standard deviation:':<35}{self.concordance.std:>10}\n"\
               f"{'Concordance standard error:':<35}{self.concordance.sem:>10}\n"

    def __repr__(self):
        return f"<DiversityCalculator(sname={self.sname}), "\
               f"concordance={self.concordance}>"


    ######################################
    ########### Private methods ##########
    ######################################

    def _calculate_concordance_stats(self, weights_col='Read_RPKM_wt'):
        """
        Based on self.reference_hits dataframe computes concordance statisctics:
         - concordance
         - concordance standard deviation
         - concordance standard error of the mean
        :param weights_col: column name containing weights or None.
        If None is supplied unweighted statistics are calculated.
        :return: tupple with (mean, std, sem)
        """
        array = np.array(self.reference_hits.iden)
        weights = np.array(self.reference_hits[weights_col]) if weights_col else np.ones(len(self.reference_hits))
        concrordance_stats = DescrStatsW(self.reference_hits.iden, weights, ddof=0)
        concordance_sem = concrordance_stats.std / (len(self.reference_hits)) ** 0.5

        return (concrordance_stats.mean, concrordance_stats.std, concordance_sem)

    ######################################################
    ##########Printing Methods############################
    ######################################################

    def print_reference_hits_summary(self):
        with pd.option_context('display.max_rows', None, 'display.max_columns',
                               None):  # more options can be specified also
            print(self.reference_hit_summary)

    ######################################################
    ###########Private Methods############################
    ######################################################

    def _summarize_reference_hits(self, df):

        def summarize_output(x):
            d = {}
            d['Ref_Size'] = x['Ref_size'].mean()
            d['Total_Reads'] = x['suid'].count()
            d['Total_aligned_nt'] = x['alen'].sum()
            d['Total_RPKM'] = x['Read_RPKM'].sum()
            d['Theoretical_Copies'] = d['Total_aligned_nt'] / d['Ref_Size']
            d['Total_RPKM_wt'] = x['Read_RPKM_wt'].sum()
            d['Average_iden'] = x['iden'].mean()
            d['Iden_VAR'] = 0 if d['Total_Reads'] == 1 else x['iden'].var()
            d['Iden_STD'] = d['Iden_VAR'] ** 0.5
            d['Iden_SEM'] = d['Iden_STD'] / d['Total_Reads'] ** 0.5

            return pd.Series(d, index=list(d.keys()))

        return df.groupby('suid').apply(summarize_output)

    def _convert_iden_percent_to_frac(self):
        """
        MODIFIES: self.reference_hits attribute:
            iden column of the reference_hits dataframe is converted to fractions from percents
        :return: void
        """
        self.reference_hits['iden'] = self.reference_hits['iden'] / 100.0

    @staticmethod
    def _load_hits(file_path):
        """
        Reads hit summary from blast runner
        :param fnames_lst: list of input filenames
        :return: dataframe with combined data from fnames_lst or None if no hits were found
        """
        df_lst = []
        empty_lst = []

        hit_data = pd.read_csv(file_path, sep='\t')
        return hit_data.drop("Unnamed: 0", axis=1)

    def calculate_novelty(self):
        return len(self.reference_hits[self.reference_hits.iden < self.concordance.mean])/len(self.reference_hits)

