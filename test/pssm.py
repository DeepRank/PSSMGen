from pssmgen.pssm import PSSMdecoy


gen = PSSMdecoy('1AK4')

# generates the FASTA query
gen.get_fasta()

# generates the PSSM
#gen.get_pssm()

# map the pssm to the pdb
gen.map_pssm()
