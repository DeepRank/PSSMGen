from pssmgen.pssm import PSSM

gen = PSSM('1AK4')

# set blast excuete and database
gen.configure(blast_exe='/home/software/blast/bin/psiblast',
            database='/data/DBs/blast_dbs/nr_v20180204/nr')

# generates the FASTA query
 gen.get_fasta()

# generates the PSSM
 gen.get_pssm(run=False)

# map the pssm to the pdb
gen.map_pssm()
