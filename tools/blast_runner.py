"""
Placeholder for the blast runner
It'll rely heavily on gnu parallel version, as it shows nice results
also it's a good idea to use .ncbirc for the database as in here
https://gtpb.github.io/PPB18/assets/15_Running-BLAST_sys.argv
:
- parallel blasting independent of tab file processing
!!! Note processing could be moved to a separate file
- processing tab-files separately to avoid large memory overheads (maybe during blasting if possible, ok sequentially)
"""

#p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')

# if os.workdir_path.exists("0.tab"):
#     os.remove("0.tab")
# cmd = "hmmscan " + "-o /dev/null --noali --cpu 1 -E " + str(
#     e_val) + " --tblout >(tee -a 0.tab)  " + model + " " + fnames
# args = shlex.split(cmd)
# # os.system(cmd)
# p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
# # print(p.stdout.read())
#
# output = p.stdout.read()
import os
from shlex import split
from subprocess import Popen, PIPE

from config import E_VAL
from tools.ext_process import run_external_command


class BlastRunner():
    def __init__(self, reference_path, workdir_path,
                 blast_dir_name="blast",
                 blast_db_dir_name=".blast_db", blast_db_name="db", e_val=E_VAL,
                 column_headers="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
        self.blast_dir_name = blast_dir_name
        self.blast_db_name = blast_db_name
        self.blast_db_dir_name = blast_db_dir_name
        self.blast_dir_path = os.path.join(workdir_path, self.blast_dir_name)
        self.reference_path = reference_path
        self.db_path = self._init_blast_dir(self.blast_dir_path, self.blast_db_dir_name)
        self.e_val = str(e_val)
        self.db_path_w_name = None #is set in _make_blast_db() if succesful
        self.column_headers = column_headers
        self._make_blast_db(self.reference_path, self.db_path)


    def run_blast_parallel(self, input_path, sname):
        out_path = os.path.join(self.blast_dir_path, sname+".tab")
        cmd = f"cat {input_path} | parallel --block 10M --recstart '>' --pipe \
        blastn  -outfmt \\'\"6 {self.column_headers} \\'\" -evalue {self.e_val} -db {self.db_path_w_name}  -query -  > {out_path}"
        run_external_command(cmd)
        return out_path

    def _init_blast_dir(self, blast_dir_path, blast_db_dir_name):
        blast_db_dir = os.path.join(blast_dir_path, blast_db_dir_name)

        if os.path.exists(blast_dir_path):
            print("Blast directory exists")
            # shutil.rmtree(blast_db_Dir)
        else:
            try:
                os.mkdir(blast_dir_path)
                os.mkdir(blast_db_dir)
            except OSError:
                print("ERROR: Cannot create project directory: " + blast_db_dir)
                raise

        return blast_db_dir

    def _make_blast_db(self, reference_path, db_path):
        db_path = os.path.join(db_path, self.blast_db_name)
        build_blast_db_cmd = ['makeblastdb', '-in', reference_path, '-dbtype', 'nucl', '-out', db_path , '-parse_seqids']
        process = Popen(build_blast_db_cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise RuntimeError(f"make_blast_db() failed with code {process.returncode}:\n {stdout} {stderr}")
        else:
            self.db_path_w_name = db_path
            return 0

