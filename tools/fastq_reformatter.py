import os

from tools.ext_process import run_external_command
from tools.fs import init_dir
from tools.utils import remove_suffix, remove_common_substring


class FastqReformatter():
    """
    Converts Fastq to fasta
    Combines sample files into one
    renames read names
    """
    # TODO optimize parser to make read renaming unnecessary
    def __init__(self, workdir_path, out_fasta_dir_name=".rf_samples", gzipped=False):
        self.out_fasta_dir_name = out_fasta_dir_name
        self.gzipped = gzipped
        self.out_fasta_dir_path = init_dir(workdir_path, self.out_fasta_dir_name, force_del=True)

    def run_fastq_reformatter(self, input_path_lst, sname):

        if self.gzipped:
            out_path = os.path.join(self.out_fasta_dir_path, sname + ".rf.fasta.gz")
        else:
            out_path = os.path.join(self.out_fasta_dir_path, sname + ".rf.fasta")
        #reformat.sh in=$mg out=stdout.fasta | rename.sh in=stdin.fasta out=$name.fasta prefix=$name

        for input_path in input_path_lst:
            filename = os.path.basename(input_path)
            prefix = remove_common_substring(filename)
            cmd = f"reformat.sh -Xmx2g in={input_path} out=stdout.fasta | rename.sh app=t in=stdin.fasta out=stdout.fasta prefix={prefix} out={out_path}"
            if self.gzipped:
                cmd = f"reformat.sh -Xmx2g in={input_path} out=stdout.fasta | rename.sh app=t in=stdin.fasta out=stdout.fasta prefix={prefix} out={out_path}"

            run_external_command(cmd)

        return out_path

