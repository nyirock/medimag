import os
from glob import glob
import csv
import pandas as pd

from config import BLAST_DIR_NAME, BLAST_DB_NAME

#TODO: use single function once the result congruence is verified
def write_df_to_csv(df, path):
    df.to_csv(path, sep='\t')

def load_blast_hits_to_df(filepath):
    blast_cols = ['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits']
    try:
        df = pd.read_csv(filepath, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=blast_cols)  # making an empty dataframe
    df.columns = blast_cols
    return df


def init_blast_dir(path, blast_dir_name, blast_db_dir_name='blast_db'):
    blast_dir = os.path.join(path, blast_dir_name)
    blast_db_dir = os.path.join(blast_dir, blast_db_dir_name)

    if os.path.exists(blast_dir):
            print("Blast directory exists")
            #shutil.rmtree(blast_db_Dir)
    else:
        try:
            os.mkdir(blast_dir)
            os.mkdir(blast_db_dir)
        except OSError:
            print("ERROR: Cannot create project directory: " + blast_db_dir)
            raise

    return blast_db_dir


def get_file_list(path, ext_mask):
    files = []
    for file in glob(os.path.join(path, ext_mask)):
        files.append(os.path.basename(file))
    return files


def write_sample_file_to_disk(df, path):
    df.to_csv(os.path.join(path, "sample_file.tsv"), sep='\t', index=False)


def load_sample_file_from_disk(path):
    return pd.read_csv(os.path.join(path, "sample_file.tsv"), sep='\t')


def load_sample_from_disk_to_dict(path):
    d = dict()
    with open(os.path.join(path, 'sample_file.tsv'), mode='r') as infile:
        reader = csv.reader(infile.read().splitlines(), delimiter='\t')
        next(reader, None)
        for row in reader:
            d[row[0]] = row[1].split(",")
    return d
