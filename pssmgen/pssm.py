import os
from pdb2sql.pdb2sqlcore import pdb2sql
import shutil
from Bio.Blast.Applications import NcbipsiblastCommandline
from .map_pssm2pdb import write_mapped_pssm_pdb


class PSSMdecoy(object):

    def __init__(self,caseID,pdbdir='pdb'):

        self.caseID = caseID
        self.pdbdir = pdbdir
        self.pdbs = os.listdir(os.path.join(self.caseID,self.pdbdir))

        self.One2ThreeDict = {
        'A' : 'ALA', 'R' : 'ARG', 'N' : 'ASN', 'D' : 'ASP', 'C' : 'CYS', 'E' : 'GLU', 'Q' : 'GLN',
        'G' : 'GLY', 'H' : 'HIS', 'I' : 'ILE', 'L' : 'LEU', 'K' : 'LYS', 'M' : 'MET', 'F' : 'PHE',
        'P' : 'PRO', 'S' : 'SER', 'T' : 'THR', 'W' : 'TRP', 'Y' : 'TYR', 'V' : 'VAL',
        'B' : 'ASX', 'U' : 'SEC', 'Z' : 'GLX'
        }

        self.Three2OneDict = {v: k for k, v in self.One2ThreeDict.items()}


    def get_fasta(self,chain=['A','B'],outdir='./fasta/'):

        outdir = os.path.join(self.caseID,outdir)
        if not os.path.isfile(outdir):
            os.mkdir(outdir)

        pdb = os.path.join(os.path.join(self.caseID,self.pdbdir),self.pdbs[0])
        sqldb = pdb2sql(pdb)
        for c in chain:

            # get the unique residues
            res = sqldb.get_residues(chainID=c)

            # get the one letter resiude
            seq = ''
            for r in res:
                seq += self.Three2OneDict[r[1]]

            # write the file
            fname = os.path.join(outdir,self.caseID + '_%s' %c + '.fasta')
            f = open(fname,'w')
            f.write('>%s' %self.caseID + '_%s' %c)
            f.write(seq)
            f.close()
        sqldb.close()


    def get_pssm(self,fasta_dir='./fasta/',
                 blast = '/home/clgeng/software/blast/bin/psiblast',
                 db = '/data/lixue/DBs/blast_dbs/nr_v20180204/nr',
                 outdir='./pssm_raw/'):

        fasta_dir = os.path.join(self.caseID,fasta_dir)
        outdir = os.path.join(self.caseID,outdir)
        if not os.path.isfile(outdir):
            os.mkdir(outdir)

        for q in os.listdir(fasta_dir):

            query = os.path.join(fasta_dir,q)
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

            name = os.path.split(os.path.basename(query))[0] + '.pssm'
            out = os.path.join(outdir,name)
            shutil.copy2(name+ '.pssm.3', out)


    def map_pssm(self,pssm_dir='pssm_raw',outdir='pssm',chain=['A','B']):

        pssm_dir = os.path.join(self.caseID,pssm_dir)
        outdir = os.path.join(self.caseID,outdir)
        if not os.path.isfile(outdir):
            os.mkdir(outdir)

        pf = os.listdir(pssm_dir)
        pssm_files = {}
        for c in chain:
            pssm_files[c] = list(filter(lambda x: x.endswith(c+'.pssm'),pf))[0]

        for p in self.pdbs:
            pdb = os.path.join(os.path.join(self.caseID,self.pdbdir),p)
            for c in chain:
                pssm = os.path.join(pssm_dir,pssm_files[c])
                print(pssm)
                write_mapped_pssm_pdb(pssm, pdb, c, outdir)

