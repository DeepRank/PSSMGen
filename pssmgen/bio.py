from Bio.Blast.Applications import NcbipsiblastCommandline
exec = '/home/clgeng/software/blast/bin/psiblast'
db = ' /data/lixue/DBs/blast_dbs/nr_v20180204/nr'
query = '1AK4/fasta/1AK4_A.fasta'
psi_cline = NcbipsiblastCommandline(exec,db = db, query = query, evalue = 0.0001,
                                    out = 'test_psi.xml', outfmt = 7, out_pssm = 'test_pssm',
				    save_each_pssm=True,num_iterations=3,save_pssm_after_last_round=True,
                                   out_ascii_pssm='1AK4.A.pssm',matrix='BLOSUM62')
print(str(psi_cline))
print(psi_cline._validate())
psi_cline()
