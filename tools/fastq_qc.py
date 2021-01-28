import os

from tools.ext_process import run_external_command
from tools.fs import init_dir
from tools.util import remove_suffix


class FastqQC():
    """
    Performs Fastq quality control
    """
    def __init__(self, workdir_path, out_fasta_dir_name=".qc_samples", gzipped=False):
        self.out_fasta_dir_name = out_fasta_dir_name
        self.out_fasta_dir_path = init_dir(workdir_path, self.out_fasta_dir_name, force_del=True)
        self.gzipped = gzipped

    def run_fastq_qc(self, input_path_lst, sname):
        out_path_lst = []
        assert 0< len(input_path_lst) <= 2, "Max number of files processed together is 2"
        # if self.gzipped:
        #     out_path = os.path.join(self.out_fasta_dir_path, sname + ".rf.fasta.gz")
        # else:
        #     out_path = os.path.join(self.out_fasta_dir_path, sname + ".rf.fasta")

        for input_path in input_path_lst:
            filename = os.path.basename(input_path)
            prefix = remove_suffix(filename, ".fastq.gz")
            cmd = f"reformat.sh in={input_path} out=stdout.fasta | rename.sh in=stdin.fasta out=stdout.fasta prefix={prefix} >> {out_path}"
            run_external_command(cmd)

        return out_path_lst