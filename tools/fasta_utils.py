import os
import shutil
import tempfile

from Bio import SeqIO
import gzip
from tools.ext_process import run_external_command
from tools.fs import is_gzipped


def parse_contigs_ind(filepath):
    """
    Returns sequences index from the input files(s)
    remember to close index object after use
    """
    if filepath.endswith("gz") or is_gzipped(filepath):
        with tempfile.NamedTemporaryFile(mode='wt') as tf:
            with gzip.open(filepath, 'rt') as f_gz:
                shutil.copyfileobj(f_gz, tf)
            temp_file_path = os.path.join(tempfile.gettempdir(), tf.name)
            record_dict = SeqIO.index(temp_file_path, "fasta")
    else:
        with open(filepath):
            record_dict = SeqIO.index(filepath, "fasta")
    return record_dict

def write_recruited_reads_to_fasta(df, all_records, outfile):
    records = []
    ids = df['quid'].tolist()
    # if len(ids)==len(sequences):
    for j in range(len(ids)):
        records.append(all_records[ids[j]])
    with open(outfile, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

def get_nt_size(record_dict):
    """
    Gets the amount of nucleotides in a file.
    Receives record_dict.
    """
    return sum(len(seq) for seq in record_dict.values())

def get_nt_size_ext(file_path):
    cmd = "seqkit stats -T "+file_path
    output = run_external_command(cmd).split()
    n_seq = int(float(output[-5]))
    nt_sz = int(float(output[-4]))
    return (n_seq, nt_sz)

def retrive_sequence(contig_lst, rec_dic):
    """
    Returns list of sequence elements from dictionary/index of SeqIO objects specific to the contig_lst parameter
    """
    contig_seqs = list()
    #record_dict = rec_dic
    #handle.close()
    for contig in contig_lst:
        contig_seqs.append(str(rec_dic[contig].seq))#fixing BiopythonDeprecationWarning
    return contig_seqs
