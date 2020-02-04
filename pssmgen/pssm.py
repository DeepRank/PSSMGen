import os, glob, shutil
from Bio.Blast.Applications import NcbipsiblastCommandline
from pdb2sql import pdb2sql

from pssmgen.map_pssm2pdb import write_mapped_pssm_pdb


class PSSM(object):

    def __init__(self,caseID='',pdbdir='pdb'):

        """Compute the PSSM and map the the sequence for a series of decoys.

        Args:
            caseID (TYPE): Name of the case. This must correspond to a directory name
            pdbdir (str, optional): directory name where the decoys pdbs are stored within caseID

        Example:

        >>> # import the module
        >>> from pssmgen import PSSM
        >>>
        >>> # create ab instance
        >>> gen = PSSM('1AK4',pdbdir='water')
        >>>
        >>> # configure the generator
        >>> gen.configure(blast='/home/software/blast/bin/psiblast',
                         database='/data/DBs/blast_dbs/nr_v20180204/nr')
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

        self.psiblast_parameter = {
        0 : { 'wordSize':2, 'gapOpen':9,  'gapExtend':1, 'scoringMatrix':'PAM30' },
        1 : { 'wordSize':3, 'gapOpen':9,  'gapExtend':1, 'scoringMatrix':'PAM30' },
        2 : { 'wordSize':3, 'gapOpen':10, 'gapExtend':1, 'scoringMatrix':'PAM70' },
        3 : { 'wordSize':3, 'gapOpen':10, 'gapExtend':1, 'scoringMatrix':'BLOSUM80'},
        4 : { 'wordSize':3, 'gapOpen':11, 'gapExtend':1, 'scoringMatrix':'BLOSUM62'}
        }


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


    def configure(self, blast_exe=None, database=None,
        num_threads = 4, evalue=0.0001, comp_based_stats='T',
        max_target_seqs=2000, num_iterations=3, outfmt=7,
        save_each_pssm=True, save_pssm_after_last_round=True):
        """Configure the blast executable, database and psiblast parameters.

        Notes:
            For more details about psiblast paramters, check 'psiblast -help'.

        Args:
            blast_exe (str): Path to the psiblast executable
            database (str) : Path to the Blast database
            num_threads (int): Number of threads (CPUs) to use in the BLAST search
            evalue (float): Expectation value (E) threshold for saving hits
            comp_based_stats (str, int): Use composition-based statistics,
                    0, F or f: no composition-based statistics
                    2, T or t, D or d : Composition-based score adjustment
                    as in Bioinformatics 21:902-911, 2005, conditioned on
                    sequence properties
            max_target_seqs (int): Maximum number of aligned sequences to keep,
                    not applicable for outfmt <= 4
            num_iterations (int): Number of iterations to perform,
                    0 means run until convergence
            save_each_pssm (bool): Save PSSM after each iteration
            save_pssm_after_last_round (bool): Save PSSM after the last database search
        """

        self.blast_exe = blast_exe
        self.blast_config = {
            'db': database,
            'num_threads': num_threads,
            'evalue': evalue,
            'comp_based_stats': comp_based_stats,
            'max_target_seqs': max_target_seqs,
            'num_iterations': num_iterations,
            'save_each_pssm': save_each_pssm,
            'save_pssm_after_last_round': save_pssm_after_last_round,
        }

    def get_pssm(self,fasta_dir='fasta/',
                 outdir='pssm_raw/',
                 run=True):

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

        out_fmt = '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi   ,\
                   sacc saccver sallacc slen qstart qend sstart send qseq sseq  ,\
                   evalue bitscore score length pident nident mismatch positive ,\
                   gapopen gaps ppos frames qframe sframe btop staxids stitle   ,\
                   salltitles sstrand qcovs qcovhsp qcovus'

        for q in os.listdir(fasta_dir):

            # get the fasta quey
            query = os.path.join(fasta_dir,q)
            name = os.path.splitext(os.path.basename(query))[0]

            # set up the output names
            out_ascii_pssm = os.path.join(outdir,name + '.pssm')
            out_pssm = os.path.join(outdir,name + '.cptpssm')
            out_xml = os.path.join(outdir,name + '.xml')

            # get the parameters
            blast_param = self._get_psiblast_parameters(query)

            # set up the psiblast calculation
            psi_cline = NcbipsiblastCommandline(
                               cmd = self.blast_exe,
                               query = query,
                               word_size = blast_param['wordSize'],
                               gapopen = blast_param['gapOpen'],
                               gapextend = blast_param['gapExtend'],
                               matrix = blast_param['scoringMatrix'],
                               outfmt = out_fmt,
                               out_ascii_pssm = out_ascii_pssm,
                               out_pssm = out_pssm,
                               out = out_xml,
                               **self.blast_config
                               )

            # check that it's correct
            psi_cline._validate()

            if run:

                # run the blast query
                psi_cline()

                # copyt the final pssm to its final name
                shutil.copy2(out_ascii_pssm + '.%d' %self.blast_config['num_iterations'], out_ascii_pssm)

                # remove all the other files
                for filename in glob.glob(out_pssm+'.*'):
                    os.remove(filename)
                for filename in glob.glob(out_ascii_pssm+'.*'):
                    os.remove(filename)
                os.remove(out_xml)


    def _get_psiblast_parameters(self,fasta_query):

        f = open(fasta_query)
        data =f.readlines()
        f.close()

        seq = 0
        for l in data[1:]:
            seq += len(l)

        if seq < 30:
            return self.psiblast_parameter[0]
        elif seq < 35:
            return self.psiblast_parameter[1]
        elif seq < 50:
            return self.psiblast_parameter[2]
        elif seq < 85:
            return self.psiblast_parameter[3]
        else:
            return self.psiblast_parameter[4]

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
                write_mapped_pssm_pdb(pssm, pdb, c, outdir)

        # write mapped pdb
        self._write_mapped_pdb(pdbindir=outdir.split('/')[1])

    def _write_mapped_pdb(self, pdbindir='pssm', pdbrawdir='pdb_raw'):
        """Write mapped pdb to working folder, by default 'pdb'

        Args:
            pdbindir (str): directory name where mapped pdb and pssm are stored
            pdbrawdir (str): directory where backup original pdb files
        """

        pdbindir = os.path.join(self.caseID, pdbindir)
        pdbrawdir = os.path.join(self.caseID, pdbrawdir)
        print('='*80)
        print(self.caseID)
        print(pdbindir, pdbrawdir)

        pdb_files = [f for f in os.listdir(pdbindir) if f.endswith('.pdb')]

        if pdb_files:
            if not os.path.isdir(pdbrawdir):
                os.mkdir(pdbrawdir)

            pdb_dict = {}
            for pdb in pdb_files:
                id, chain, _, _ = pdb.split('.')
                if id not in pdb_dict:
                    pdb_dict[id] = []
                pdb_dict[id].append(chain)

            for id, chains in pdb_dict.items():
                pdb_wd = os.path.join(self.caseID, self.pdbdir, id+".pdb")
                pdb_raw = os.path.join(pdbrawdir, id+".pdb")
                os.rename(pdb_wd, pdb_raw)
                if len(chains) == 1:
                    pdb_new = os.path.join(pdbindir, id + "." + chains[0] + '.pssm.pdb')
                    os.rename(pdb_new, pdb_wd)
                else:
                    with open(pdb_wd, 'w') as f:
                        # write REMARK
                        pdb_new = os.path.join(pdbindir, id + "." + chains[0] + '.pssm.pdb')
                        with open(pdb_new, 'r') as fpdb:
                            lines = [line for line in fpdb if line.startswith('REMARK')]
                            f.writelines(lines)

                        # write each chain
                        for chain in chains:
                            pdb_new = os.path.join(pdbindir, id + "." + chain + '.pssm.pdb')
                            with open(pdb_new, 'r') as fpdb:
                                lines = [line for line in fpdb
                                    if line.startswith('ATOM') and
                                    line[21:22] == chain]
                                f.writelines(lines)
                            os.remove(pdb_new)
                        f.write('END\n')
