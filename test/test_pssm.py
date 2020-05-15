from pssmgen.pssm import PSSM

# initiate the PSSM object
gen = PSSM(work_dir='7CEI')

# set psiblast executable, database and other psiblast parameters (here shows the defaults)
gen.configure(blast_exe='/home/software/blast/bin/psiblast',
            database='/data/DBs/blast_dbs/nr_v20180204/nr',
            num_threads = 4, evalue=0.0001, comp_based_stats='T',
            max_target_seqs=2000, num_iterations=3, outfmt=7,
            save_each_pssm=True, save_pssm_after_last_round=True)

# generates FASTA files
gen.get_fasta(pdb_dir='pdb', chain=('A','B'), out_dir='fasta')

# generates PSSM
gen.get_pssm(fasta_dir='fasta', out_dir='pssm_raw', run=True)

# map PSSM and PDB to get consisitent files
gen.map_pssm(pssm_dir='pssm_raw', pdb_dir='pdb', out_dir='pssm', chain=('A','B'))

# write consistent files and move
gen.get_mapped_pdb(pdbpssm_dir='pssm', pdb_dir='pdb', pdbnonmatch_dir='pdb_nonmatch')