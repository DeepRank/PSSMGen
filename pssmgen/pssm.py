import os
import glob
import shutil
import re
import logging
from Bio.Blast.Applications import NcbipsiblastCommandline
from pdb2sql import pdb2sql

from pssmgen.map_pssm2pdb import write_mapped_pssm_pdb

logger = logging.getLogger()

class PSSM():

    def __init__(self, work_dir='.'):
        """Compute PSSM and map the sequence for a series of decoys.

        Args:
            work_dir (str, optional): the working directory that contains the
                pdb directory and/or fasta directory.
                Defaults to '.'.
        """

        self.work_dir = work_dir

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


    def get_fasta(self, pdb_dir='pdb', pdbref=None, chain=('A','B'), out_dir='fasta'):
        """Extract the sequence of the chains and writes a fasta query file for each.

        Args:
            pdb_dir (str, optional): path for pdb files.
                Defaults to 'pdb'.
            pdbref (str, optional): filename of the reference pdb.
                Defaults to None, i.e. randomly select one.
            chain (tuple, optional): Name of the chains in the pdbs.
                Defaults to ('A','B').
            out_dir (str, optional): output path for fasta files.
                Defaults to 'fasta'.
        """
        out_dir = os.path.join(self.work_dir,out_dir)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        if pdbref:
            pdb = os.path.join(self.work_dir, pdb_dir, pdbref)
        else:
            pdbs = os.listdir(os.path.join(self.work_dir, pdb_dir))
            pdbs = list(filter(lambda x: x.endswith('.pdb'), pdbs))
            pdb = os.path.join(self.work_dir, pdb_dir, pdbs[0])

        sqldb = pdb2sql(pdb)
        print('Generating FASTA files...')
        logger.info('Output FASTA files:')
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
            caseID = re.split('_|\.', os.path.basename(pdb))[0]
            fname = os.path.join(out_dir, caseID + '.%s' %c + '.fasta')
            with open(fname,'w') as f:
                f.write('>%s' %caseID + '.%s\n' %c)
                f.write(seq)
            logger.info(f'  {fname}')

        print(f'FASTA files generated in {out_dir}.\n')

    def configure(self, blast_exe=None, database=None,
        num_threads = 4, evalue=0.0001, comp_based_stats='T',
        max_target_seqs=2000, num_iterations=3, outfmt=7, outspecifiers = None,
        save_each_pssm=True, save_pssm_after_last_round=True):
        """Configure the blast executable, database and psiblast parameters.

        Notes:
            For more details about psiblast paramters, check 'psiblast -help'.

        Args:
            blast_exe (str): Path to the psiblast executable
            database (str) : Path to the Blast database
            num_threads (int): Number of threads (CPUs) to use in the BLAST search.
                Defaults to 4.
            evalue (float): Expectation value (E) threshold for saving hits.
                Defaults to 0.0001.
            comp_based_stats (str, int): Use composition-based statistics,
                    0, F or f: no composition-based statistics
                    2, T or t, D or d : Composition-based score adjustment
                    as in Bioinformatics 21:902-911, 2005, conditioned on
                    sequence properties
                Defaults to 'T'.
            max_target_seqs (int): Maximum number of aligned sequences to keep,
                    not applicable for outfmt <= 4.
                    Defaults to 2000.
            num_iterations (int): Number of iterations to perform,
                    0 means run until convergence.
                    Defaults to 3.
            outfmt (int): Format for output alignment view. Valid values are 0-18.
                    Default format 7 is "tabular with comment lines".
                    Check `psiblast -help` for details.
            outspecifiers (str): Output format 6, 7 and 10 can be additionally
                    configured to produce a custom format specified by space
                    delimited format specifiers. By default using all sepcifiers.
                    Check `psiblast -help` for details.
            save_each_pssm (bool): Save PSSM after each iteration.
                    Defaults to True.
            save_pssm_after_last_round (bool): Save PSSM after the last database search.
                    Defaults to True.
        """
        # output format sepcifiers for format 6, 7 and 10. See psiblast for details.
        specifiers = 'qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi \
                sacc saccver sallacc slen qstart qend sstart send qseq sseq \
                evalue bitscore score length pident nident mismatch positive \
                gapopen gaps ppos frames qframe sframe btop staxids stitle \
                salltitles sstrand qcovs qcovhsp qcovus'

        outfmt = str(outfmt)
        if outfmt in ['6', '7', '10']:
            if outspecifiers is not None:
                outfmt = outfmt + ' ' + outspecifiers
            else:
                outfmt = outfmt + ' ' + specifiers

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
            'outfmt': outfmt,
        }

    def get_pssm(self, fasta_dir='fasta',
                 out_dir='pssm_raw',
                 run=True,
                 save_all_psiblast_output=False):
        """Compute the PSSMs by running BLAST

        Args:
            fasta_dir (str, optional): path of fasta files.
                Defaults to 'fasta'.
            out_dir (str, optional): output path for pssm files.
                Defaults to 'pssm_raw'.
            run (bool, optional): run blast.
                Defaults to True.
            save_all_psiblast_output (bool, optional):
                save all output from psiblast, including pssm files of
                each round and details of homologs, etc.
                Defaults to False.

        Raises:
            FileNotFoundError: BLAST failed to find homologs.
        """
        fasta_dir = os.path.join(self.work_dir,fasta_dir)
        out_dir = os.path.join(self.work_dir,out_dir)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        for q in os.listdir(fasta_dir):

            # get the fasta quey
            query = os.path.join(fasta_dir,q)
            name = os.path.splitext(os.path.basename(query))[0]

            # set up the output names
            out_ascii_pssm = os.path.join(out_dir,name + '.pssm')
            out_pssm = os.path.join(out_dir,name + '.cptpssm')
            out_homologs = os.path.join(out_dir,name + '.homologs')

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
                               out_ascii_pssm = out_ascii_pssm,
                               out_pssm = out_pssm,
                               out = out_homologs,
                               **self.blast_config
                               )

            # check that it's correct
            psi_cline._validate()

            if run:
                # run the blast query
                print('Generatinng raw PSSMs with BLAST...')
                psi_cline()

                # copy the pssm of last exiting iteration to its final name
                logger.info('Generated raw PSSM files:')
                for i in reversed(range(self.blast_config['num_iterations'])):
                    fpssm = out_ascii_pssm + '.' + str(i+1)
                    if os.path.isfile(fpssm):
                        shutil.copy2(fpssm, out_ascii_pssm)
                        logger.info(f'  {out_ascii_pssm}')
                        break
                    elif i==0:
                        raise FileNotFoundError(f'Not found {fpssm}. \
                            PSIBlast may fail to find homologs for given fasta')

                if not save_all_psiblast_output:
                    # remove all the other files
                    for filename in glob.glob(out_pssm+'.*'):
                        os.remove(filename)
                    for filename in glob.glob(out_ascii_pssm+'.*'):
                        os.remove(filename)
                    os.remove(out_homologs)
                else:
                    logger.info(f'Other intermediate output files are in {out_dir}')
                print(f'Raw PSSM files generated in {out_dir}.\n')


    def _get_psiblast_parameters(self, fasta_query):

        f = open(fasta_query)
        data =f.readlines()
        f.close()

        seq = 0
        for l in data[1:]:
            seq += len(l)

        if seq < 30:
            p = self.psiblast_parameter[0]
        elif seq < 35:
            p = self.psiblast_parameter[1]
        elif seq < 50:
            p = self.psiblast_parameter[2]
        elif seq < 85:
            p = self.psiblast_parameter[3]
        else:
            p = self.psiblast_parameter[4]

        return p

    def map_pssm(self, pssm_dir='pssm_raw', pdb_dir='pdb', out_dir='pssm', chain=('A','B')):
        """Map the raw pssm files to the pdb files of the decoys

        Args:
            pssm_dir (str, optional): path of raw pssm files.
                Defaults to 'pssm_raw'.
            pdb_dir (str, optional): path of pdb files.
                Defaults to 'pdb'.
            out_dir (str, optional): output path for mapped pssm files.
                Defaults to 'pssm'.
            chain (tuple, optional): chain names. Defaults to ('A','B').
        """
        pssm_dir = os.path.join(self.work_dir,pssm_dir)
        out_dir = os.path.join(self.work_dir,out_dir)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        # get list of pssm files
        pf = os.listdir(pssm_dir)
        pssm_files = {}
        for c in chain:
            pssm_files[c] = list(filter(lambda x: x.endswith(c+'.pssm'),pf))[0]

        # get list of pdb files
        pdbs = os.listdir(os.path.join(self.work_dir, pdb_dir))
        pdbs = list(filter(lambda x: x.endswith('.pdb'), pdbs))

        # map pssm and pdb
        print('Generatinng mapped PSSMs...')
        logger.info('\nOutput mapped PSSM files:')
        for p in pdbs:
            pdb = os.path.join(os.path.join(self.work_dir, pdb_dir), p)
            for c in chain:
                pssm = os.path.join(pssm_dir,pssm_files[c])
                write_mapped_pssm_pdb(pssm, pdb, c, out_dir)
        print(f'Mapped PSSM files generated in {out_dir}.\n')

    def get_mapped_pdb(self, pdbpssm_dir='pssm', pdb_dir='pdb',
        pdbnonmatch_dir='pdb_nonmatch'):
        """Write mapped pdb to working folder, by default 'pdb'.

        Args:
            pdbpssm_dir (str, optional): path of mapped pssm files.
                Defaults to 'pssm'.
            pdb_dir (str, optional): output path of mapped pdb files.
                Defaults to 'pdb'.
            pdbnonmatch_dir (str, optional): output path of nonconsistent pdb files.
                Defaults to 'pdb_nonmatch'.
        """
        pdbpssm_dir = os.path.join(self.work_dir, pdbpssm_dir)
        pdbnonmatch_dir = os.path.join(self.work_dir, pdbnonmatch_dir)

        pdb_files = [f for f in os.listdir(pdbpssm_dir) if f.endswith('.pdb')]

        if pdb_files:
            print('Generatinng mapped PDB files...')
            if not os.path.isdir(pdbnonmatch_dir):
                os.mkdir(pdbnonmatch_dir)

            pdb_dict = {}
            for pdb in pdb_files:
                caseid, chain, _, _ = pdb.split('.')
                if caseid not in pdb_dict:
                    pdb_dict[caseid] = []
                pdb_dict[caseid].append(chain)

            for caseid, chains in pdb_dict.items():
                pdb_wd = os.path.join(self.work_dir, pdb_dir, caseid+".pdb")
                pdb_raw = os.path.join(pdbnonmatch_dir, caseid+".pdb")
                os.rename(pdb_wd, pdb_raw)
                if len(chains) == 1:
                    pdb_new = os.path.join(pdbpssm_dir, caseid + "." + chains[0] + '.pssm.pdb')
                    os.rename(pdb_new, pdb_wd)
                else:
                    with open(pdb_wd, 'w') as f:
                        # write REMARK
                        pdb_new = os.path.join(pdbpssm_dir, caseid + "." + chains[0] + '.pssm.pdb')
                        with open(pdb_new, 'r') as fpdb:
                            lines = [line for line in fpdb if line.startswith('REMARK')]
                            f.writelines(lines)

                        # write each chain
                        for chain in chains:
                            pdb_new = os.path.join(pdbpssm_dir, caseid + "." + chain + '.pssm.pdb')
                            with open(pdb_new, 'r') as fpdb:
                                lines = [line for line in fpdb
                                    if line.startswith('ATOM') and
                                    line[21:22] == chain]
                                f.writelines(lines)
                            os.remove(pdb_new)
                        f.write('END\n')
            print(f'Inconsistent raw PDB files moved to {pdbnonmatch_dir}')
            print(f'Mapped PDB files generated in {os.path.join(self.work_dir, pdb_dir)}.\n')
