import os

from tools.fs import init_dir
from tools.util import remove_common_substring
from tools.ext_process import run_external_command

class FastqQC():
    """
        Performs Fastq quality control
        """
    def __init__(self, workdir_path, out_fasta_dir_name=".qc_samples", gzipped=True):
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

        if len(input_path_lst) == 2:
            if "r1" in input_path_lst[0].lower():
                infile1, infile2 = input_path_lst[0], input_path_lst[1]
            else:
                infile1, infile2 = input_path_lst[1], input_path_lst[0]

            outfile1 = os.path.join(self.out_fasta_dir_path, remove_common_substring(os.path.basename(infile1)) + ".qc.fastq")
            outfile2 = os.path.join(self.out_fasta_dir_path, remove_common_substring(os.path.basename(infile2)) + ".qc.fastq")
            if self.gzipped:
                outfile1 = outfile1 + ".gz"
                outfile2 = outfile2 + ".gz"
            cmd = f"bbduk.sh  -Xmx2g overwrite=true in1={infile1} in2={infile2} out1={outfile1} out2={outfile2} ref=adapters.fa \
ktrim=r k=23 mink=11 hdist=2 minlen=25 qtrim=rl trimq=4 hdist=1 tpe tbo"

        elif len(input_path_lst) == 1:
            #single_end_workflow
            infile = input_path_lst[0]
            outfile = os.path.join(self.out_fasta_dir_path, remove_common_substring(os.path.basename(infile)) + ".qc.fastq")
            if self.gzipped:
                outfile = outfile + ".gz"
            cmd = f"bbduk.sh  -Xmx2g overwrite=true in={infile}  out={outfile} ref=adapters.fa \
ktrim=r k=23 mink=11 hdist=2 minlen=25 qtrim=rl trimq=4 hdist=1 tpe tbo"

        else:
            raise RuntimeError("QC expects 1 or 2 fastq files per sample")

        qc_report = run_external_command(cmd)
        with open(os.path.join(self.out_fasta_dir_path, "bbduk_report.log"), "a") as fh:
            fh.write(qc_report.decode())

        if len(input_path_lst) == 2:
            return [outfile1, outfile2]
        elif len(input_path_lst) == 1:
            return [outfile]
        else:
            raise RuntimeError("QC expects 1 or 2 fastq files per sample")




