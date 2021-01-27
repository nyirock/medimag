import os
import shutil
from glob import glob
import csv
import pandas as pd

from config import BLAST_DIR_NAME, BLAST_DB_NAME

#TODO: use single function once the result congruence is verified
def write_df_to_csv(df, path):
    df.to_csv(path, sep='\t')

def write_dct_to_csv(dct, dir_path, filename, header=None):
    file_path = os.path.join(dir_path, filename)

    with open(file_path, 'w') as f:
        for key,val in dct.items():
            f.write(f"{key}\t{val}\n")

def load_blast_hits_to_df(filepath):
    blast_cols = ['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits']
    try:
        df = pd.read_csv(filepath, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=blast_cols)  # making an empty dataframe
    df.columns = blast_cols
    return df


def init_dir(path, inner_dir_name="", force_del=False):
    inner_dir_path = os.path.join(os.path.abspath(path), inner_dir_name)

    #TODO: change print statements to logging
    if os.path.exists(inner_dir_path):
            print(f"{inner_dir_name} directory exists")
            if force_del:
                shutil.rmtree(inner_dir_path)
                print(f"{inner_dir_name} directory deleted")
                mkdir(inner_dir_path)
    else:
        mkdir(inner_dir_path)
    return inner_dir_path


def mkdir(inner_dir_path):
    try:
        os.mkdir(inner_dir_path)
    except OSError:
        print("ERROR: Cannot create project directory: " + inner_dir_path)
        raise


def get_file_list(dir_path, ext_mask):
    files = []
    for file in glob(os.path.join(dir_path, ext_mask)):
        files.append(os.path.basename(file))
    return files

def write_sample_file_dct_to_disk(dct, dir_path,
                                  file_name="sample_file.tsv",
                                  header=('SampleName', 'Files')):
    file_path = os.path.join(dir_path, file_name)
    if header:
        with open(file_path, "w") as f:
            f.write(str(header[0])+"\t"+str(header[1])+"\n")

    with open(file_path, 'a') as f:
        for key, val in dct.items():
            comma_sep_vals = ','.join(val)
            f.write(f"{key}\t{comma_sep_vals}\n")

def write_sample_file_to_disk(df, path):
    df.to_csv(os.path.join(path, "sample_file.tsv"), sep='\t', index=False)


def load_sample_file_from_disk(path):
    return pd.read_csv(os.path.join(path, "sample_file.tsv"), sep='\t')


def load_sample_from_disk_to_dict(dir_path, filename="sample_file.tsv"):
    d = dict()
    with open(os.path.join(dir_path, filename), mode='r') as infile:
        reader = csv.reader(infile.read().splitlines(), delimiter='\t')
        next(reader, None)
        for row in reader:
            d[row[0]] = row[1].split(",")
    return d
