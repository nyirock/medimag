import os

BLAST_DIR_NAME = "blast"
BLAST_DB_NAME = "blast_db/db"
ALEN_BP = None
IDEN = 45
E_VAL = 1E-4
ALEN_PERCENT = 70
WORKDIR    = os.path.abspath("out")
INPUT_DIR  = os.path.abspath("in")
REFERENCE = os.path.abspath("../all_pmoA_nr99.fasta")