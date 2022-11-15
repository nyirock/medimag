"""
produce a mapping for the sample name and the raw data files
Sample names are extracted from fastq data
"""
import re
from collections import defaultdict
import pandas as pd

from tools.fs import get_file_list

def get_sample_file_dct(path, ext="*", identify_pairs=True, strict_pairs=False):
    filenames = get_file_list(path, ext_mask=ext)
    return extract_sample_names_from_filenames(filenames, identify_pairs, strict_pairs=strict_pairs)

def make_sample_file(path, ext="*"):

    filenames = get_file_list(path, ext_mask=ext)
    snames2filenames = extract_sample_names_from_filenames(filenames)

    return snames2filenames_dict_to_df(snames2filenames)


def snames2filenames_dict_to_df(dct):
    joined_filenames = {key: ','.join(val) for key, val in dct.items()}
    sample_file = pd.DataFrame(list(joined_filenames.items()), columns=["SampleName", "Files"])
    return sample_file


def extract_sample_names_from_filenames(filenames_lst, identify_pairs, strict_pairs=False):
    """
    """
    ext_mask = ['fa', 'fna', 'fq', 'fasta', 'fastq', 'fastq', 'gz', 'ac', 'qc', 'qt', 'csv', 'tab', 'rf', 'genome', 'stats','tsv']
    snames2files = defaultdict(list)
    tagged = []
    non_tagged = []

    #firts pass to identify paired reads
    if identify_pairs and not strict_pairs:
        find_paired_unbound(filenames_lst, non_tagged, snames2files, tagged)
        if ((len(non_tagged) == 0) and (len(tagged) == len(filenames_lst))):
            return snames2files
    elif identify_pairs and strict_pairs:
        find_paired_strict(filenames_lst, non_tagged, snames2files, tagged)
        if ((len(non_tagged) == 0) and (len(tagged) == len(filenames_lst))):
            return snames2files
    else:
        non_tagged = sorted(filenames_lst)

    #second pass adding single files
    find_singles(ext_mask, non_tagged, snames2files, tagged)

    assert (len(non_tagged) == 0) and (len(tagged) == len(filenames_lst)), \
        "Some input files left unprocessed"

    return dict(snames2files)


def find_singles(ext_mask, non_tagged, snames2files, tagged):
    for fname in non_tagged.copy():
        fragments = [x for x in re.split("[_.]+", fname.lower()) if x not in ext_mask]
        sname = "_".join(fragments)
        # print(sname)
        snames2files[sname].append(fname)
        tagged.append(fname)
        non_tagged.remove(fname)


def find_paired_unbound(filenames_lst, non_tagged, snames2files, tagged):
    for fname in sorted(filenames_lst):
        fragments = re.split("[_.]+r[12]?", fname.lower())
        if (len(fragments) > 1):
            lead_fragment = fragments[0]
            sname = re.split("_l00[1-9]?", lead_fragment)[0]
            tagged.append(fname)
            snames2files[sname].append(fname)
        else:
            non_tagged.append(fname)


def find_paired_strict(filenames_lst, non_tagged, snames2files, tagged):
    for fname in sorted(filenames_lst):
        fragments = re.split("[_.]+r[12]", fname.lower())
        if (len(fragments) > 1):

            sname = fragments[0]
            #sname = re.split("_l00[1-9]?", lead_fragment)[0]

            tagged.append(fname)
            snames2files[sname].append(fname)
        else:
            non_tagged.append(fname)


