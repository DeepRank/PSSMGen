# PSSMGen

Generates PSSM files for deep learning features

## Install

1. Make sure BLAST is installed and its database is available on your machine. Otherwise, install BLAST and download its databases by following the [BLAST guide](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
2. Install the PSSMgen by `pip install PSSMGen`.

## Quick Example
So far the class `PSSM` is geared toward computing the pssm files for all the decoys of a particular case. The code assumes your files have following structure:

```
 caseID
 |_ pdb
 |_ fasta
 |_ pssm_raw
 |_ pssm
```

only the `pdb` dir must exist at run time. This directory must contain the pdb files of the decoys generated for example by HADDOCK or ZDOCK. Based on these pdbs, the code will generate the fasta queries (stored in the `fasta` subdir). It will then use `psiblast` to compute the PSSM files (stored in `pssm_raw` subdir), Finally it will allign the PSSM files to the pdb sequences and store the final PSSM files in `pssm`. If the names of the subdir are different this can be specified in arguments of the method.

To use the code to generate all the pssm of caseID `1AK4` you can use:

```python
from pssmgen import PSSM

gen = PSSM('1AK4')

# set blast excuete and database
gen.configure(blast='/home/software/blast/bin/psiblast',
            database='/data/DBs/blast_dbs/nr_v20180204/nr')

# generates the FASTA query
gen.get_fasta()

# generates the PSSM
gen.get_pssm()

# map the pssm to the pdb
gen.map_pssm()
```

**If you want to use customed fasta files, you can ignore the above step `gen.get_fasta()` and mannualy create the `fasta` folder and put customed fasta files there.**

Run that on the cluster since the calculation of the PSSM files can take a few hours