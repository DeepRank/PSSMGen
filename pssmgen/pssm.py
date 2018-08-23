import os, glob, shutil
from Bio.Blast.Applications import NcbipsiblastCommandline
from pdb2sql.pdb2sqlcore import pdb2sql
from .map_pssm2pdb import write_mapped_pssm_pdb


class PSSMdecoy(object):

    def __init__(self,caseID,pdbdir='pdb'):

        """Compute the PSSM and map the the sequence for a series of decoys.

        Args:
            caseID (TYPE): Name of the case. This must correspond to a directory name
            pdbdir (str, optional): directory name where the decoys pdbs are stored within caseID

        Example:

        >>> # import the module
        >>> from pssmgen.pssm import PSSMdecoy
        >>>
        >>> # create ab instance
        >>> gen = PSSMdecoy('1AK4',pdbdir='water')
        >>>
        >>> # generates the FASTA query
        >>> gen.get_fasta()
        >>>
        >>> # generates the PSSM
        >>> gen.get_pssm()
        >>>
        >>> # map the pssm to the pdb
        >>> gen.map_pssm()

        """

        self.caseID = caseID
        self.pdbdir = pdbdir
        self.pdbs = os.listdir(os.path.join(self.caseID,self.pdbdir))
        self.pdbs = list(filter(lambda x: x.endswith('.pdb'),self.pdbs))

        self.One2ThreeDict = {
        'A' : 'ALA', 'R' : 'ARG', 'N' : 'ASN', 'D' : 'ASP', 'C' : 'CYS', 'E' : 'GLU', 'Q' : 'GLN',
        'G' : 'GLY', 'H' : 'HIS', 'I' : 'ILE', 'L' : 'LEU', 'K' : 'LYS', 'M' : 'MET', 'F' : 'PHE',
        'P' : 'PRO', 'S' : 'SER', 'T' : 'THR', 'W' : 'TRP', 'Y' : 'TYR', 'V' : 'VAL',
        'B' : 'ASX', 'U' : 'SEC', 'Z' : 'GLX'
        }

        self.Three2OneDict = {v: k for k, v in self.One2ThreeDict.items()}


    def get_fasta(self,chain=['A','B'],outdir='./fasta/'):

        """Extract the sequence of the chains and writes a fasta query file for each.

        Args:
            chain (list, optional): Name of the chains in the pdbs
            outdir (str, optional): name pf the output directory where to store the fast queries
        """

        outdir = os.path.join(self.caseID,outdir)
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        pdb = os.path.join(os.path.join(self.caseID,self.pdbdir),self.pdbs[0])
        sqldb = pdb2sql(pdb)
        for c in chain:

            # get the unique residues
            res = sqldb.get_residues(chainID=c)

            # get the one letter resiude
            seq = ''
            count = 0
            for r in res:
                seq += self.Three2OneDict[r[1]]
                count += 1
                if count == 79:
                    seq += '\n'
                    count = 0

            # write the file
            fname = os.path.join(outdir,self.caseID + '_%s' %c + '.fasta')
            f = open(fname,'w')
            f.write('>%s' %self.caseID + '_%s\n' %c)
            f.write(seq)
            f.close()
        sqldb.close()


    def get_pssm(self,fasta_dir='./fasta/',
                 blast = '/home/clgeng/software/blast/bin/psiblast',
                 db = '/data/lixue/DBs/blast_dbs/nr_v20180204/nr',
                 outdir='./pssm_raw/',
                 num_iterations=3):

        """Compute the PSSM files

        Args:
            fasta_dir (str, optional): irectory where the fasta queries are stored
            blast (str, optional): path to the psiblast executable
            db (str, optional): path to the blast database
            outdir (str, optional): output directory where to store the pssm files
            num_iterations (int, optional): number of iterations for the blast calculations
        """

        fasta_dir = os.path.join(self.caseID,fasta_dir)
        outdir = os.path.join(self.caseID,outdir)
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        for q in os.listdir(fasta_dir):

            # get the fasta quey
            query = os.path.join(fasta_dir,q)
            name = os.path.splitext(os.path.basename(query))[0]

            # set up the output names
            out_ascii_pssm = os.path.join(outdir,name + '.pssm')
            out_pssm = os.path.join(outdir,name + '.cptpssm')
            out_xml = os.path.join(outdir,name + '.xml')

            # set up the psiblast calculation
            psi_cline = NcbipsiblastCommandline(
                               blast,
                               db = db,
                               query = query,
                               evalue = 0.0001,
                               matrix='BLOSUM62',
                               outfmt = 7,
                               save_each_pssm=True,
                               num_iterations=num_iterations,
                               save_pssm_after_last_round=True,
                               out_ascii_pssm = out_ascii_pssm,
                               out_pssm = out_pssm,
                               out = out_xml
                               )
            # check that it's correct
            psi_cline._validate()

            # run the blast query
            psi_cline()

            # copyt the final pssm to its final name
            shutil.copy2(out_ascii_pssm + '.%d' %num_iterations, out_ascii_pssm)

            # remove all the other files
            for filename in glob.glob(out_pssm+'.*'):
                os.remove(filename)
            for filename in glob.glob(out_ascii_pssm+'.*'):
                os.remove(filename)
            os.remove(out_xml)

    def map_pssm(self,pssm_dir='pssm_raw',outdir='pssm',chain=['A','B']):

        """Map the raw pssm files to the pdb files of the decoys

        Args:
            pssm_dir (str, optional): name pf the directory where the pssm are stored
            outdir (str, optional): name where thmapped pssm files are stored
            chain (list, optional): name of the chains
        """
        pssm_dir = os.path.join(self.caseID,pssm_dir)
        outdir = os.path.join(self.caseID,outdir)
        if not os.path.isdir(outdir):
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

