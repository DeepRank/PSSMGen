from pdb2sql.pdb2sqlcore import pdb2sql
from .map_pssm2pdb import write_mapped_pssm_pdb
from Bio.Blast.Applications import NcbipsiblastCommandline
import shutil

class PSSM(object):

    def __init__(self,pdb,fasta,pssm):

        self.caseID = os.path.splitext(os.path.basename(pdb))[0]
        self.pdb = pdb2sql(pdb)
        self.fasta = fasta
        self.pssm = pssm

        self.One2ThreeDict = {
        'A' : 'ALA', 'R' : 'ARG', 'N' : 'ASN', 'D' : 'ASP', 'C' : 'CYS', 'E' : 'GLU', 'Q' : 'GLN',
        'G' : 'GLY', 'H' : 'HIS', 'I' : 'ILE', 'L' : 'LEU', 'K' : 'LYS', 'M' : 'MET', 'F' : 'PHE',
        'P' : 'PRO', 'S' : 'SER', 'T' : 'THR', 'W' : 'TRP', 'Y' : 'TYR', 'V' : 'VAL',
        'B' : 'ASX', 'U' : 'SEC', 'Z' : 'GLX'
        }

        self.Three2OneDict = {v: k for k, v in resmap.items()}


    def get_fasta(self,chain=['A','B']):

        for c in chain:

            # get the unique residues
            res = pdb.get_residues(chainID=c)

            # get the one letter resiude
            seq = ''
            for r in res:
                seq += self.Three2OneDict[r[1]]

            # write the file
            fname = self.caseID + '_%s' %c + '.fasta'
            f = open('fname','w')
            f.write('>%s' %self.caseID + '_%s' %c)
            f.write(seq)
            f.close()

    def get_pssm(self,blast,db,query,out):

        name = os.path.splitext(os.path.basename(query))[0]

        psi_cline = NcbipsiblastCommandline(
                           blast,
                           db = db,
                           query = query,
                           evalue = 0.0001,
                           matrix='BLOSUM62',
                           outfmt = 7,
                           save_each_pssm=True,
                           num_iterations=3,
                           save_pssm_after_last_round=True,
                           out_ascii_pssm=name + '.pssm',
                           out_pssm = name + '.cptpssm',
                           out = name + '.xml',
                           )

        shutil.copy2(name+ '.pssm.3', out)

    def map_pssm(self,fpssm,fpdb,chainID,outdir):
        write_mapped_pssm_pdb(fpssm, fpdb, chainID, outdir)

