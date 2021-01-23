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

# if os.path.exists("0.tab"):
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


class BlastRunner():
    def __init__(self, reference_path, db_path, e_val=E_VAL):
        self.reference_path = reference_path
        self.db_path = db_path
        self.e_val = str(e_val)
        self._make_blast_db(self.reference_path, self.db_path)

    def run_blast_parallel(self, input_file, out_path):
        cmd = f"cat {input_file} | parallel --block 10M --recstart '>' --pipe \
        blastn -evalue {self.e_val} -outfmt 6 -db {self.db_path}  -query - > {out_path}"
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise RuntimeError(f"run_blast_parallel() failed with code {process.returncode}:\n {stdout} {stderr}")
        else:
            #print(stdout)
            return 0

    @staticmethod
    def _make_blast_db(reference_path, db_path):
        build_blast_db_cmd = ['makeblastdb', '-in', reference_path, '-dbtype', 'nucl', '-out', db_path , '-parse_seqids']
        process = Popen(build_blast_db_cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise RuntimeError(f"make_blast_db() failed with code {process.returncode}:\n {stdout} {stderr}")
        else:
            #print(stdout)
            return 0

