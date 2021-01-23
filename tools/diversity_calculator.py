from collections import namedtuple
from typing import List
import pandas as pd
import numpy as np
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW

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
    def __init__(self, fnames_lst: List, sample_name: str):
        """
        Import csv-files as a dataframe and merge the resultant dataframe.
        :param fnames_lst:  List of csv-filenames to process
        :param sample_name: String file name
        """
        self.no_hits = False
        self.sample_name = sample_name
        self.reference_hits = self._load_hits(fnames_lst)
        if not len(self.reference_hits):
            self.no_hits = True
            #TODO: implement better the case for no hits
            return
        self._convert_iden_percent_to_frac()
        self.reference_hit_summary = self._summarize_reference_hits()
        self.concordance = Concordance(*self._calculate_concordance_stats())

    ######################################
    ########### Overloaded methods #######
    ######################################

    def __str__(self):
        return f"Sample name: {self.sample_name}\n" \
               f"\tConcordance Statistics:\n" \
               f"{'Concordance mean:':<35}{self.concordance.mean:>10}\n" \
               f"{'Concordance standard deviation:':<35}{self.concordance.std:>10}\n"\
               f"{'Concordance standard error:':<35}{self.concordance.sem:>10}\n"

    def __repr__(self):
        return f"<DiversityCalculator(sample_name={self.sample_name}), "\
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


    def _summarize_reference_hits(self):
        pass

    def _convert_iden_percent_to_frac(self):
        """
        MODIFIES: self.reference_hits attribute:
            iden column of the reference_hits dataframe is converted to fractions from percents
        :return: void
        """
        self.reference_hits['iden'] = self.reference_hits['iden'] / 100.0

    @staticmethod
    def _load_hits(fnames_lst: List):
        """
        Reads hit summary from blast runner
        :param fnames_lst: list of input filenames
        :return: dataframe with combined data from fnames_lst or None if no hits were found
        """
        df_lst = []
        empty_lst = []

        for fname in fnames_lst:
            hit_data = pd.read_csv(fname, sep='\t')
            if(len(hit_data)):
                df_lst.append(hit_data)
            else:
                empty_lst.append(hit_data)
        if df_lst:
            return pd.concat(df_lst, axis=0).drop("Unnamed: 0",axis=1)
        else:
            return pd.concat(empty_lst, axis=0).drop("Unnamed: 0", axis=1)
