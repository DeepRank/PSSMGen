# PSSMGen

Generates PSSM files for deep learning features

## Install

The only dependency is `pdb2sql`.Then clone the repo and use `pip install -e ./` as usual

## Use
So far the class `PSSMdecoy` is geared toward computing the pssm files for all the decoy of a particular case. The code assumes the following structure of the data:

 caseID
 |_ pdb
 |_ fasta
 |_ pssm_raw
 |_ pssm

only the pdb dir must exist at run time. This directory must contain the pdb files of the decoys generated for example by HADDOCK or ZDOCK. Based on these pdbs, the code will generate the fasta queries (stored in the `fasta` subdir). It will then use `psiblast` to compute the PSSM files (stored in `pssm_raw` subdir), Finally it will allign the PSSM files to the pdb sequences and store the final PSSM files in `pssm`.

To use the code to generate all the pssm of caseID `1AK4` you can use:

```python
from pssmgen.pssm import PSSMdecoy

gen = PSSMdecoy('1AK4')

# generates the FASTA query
#gen.get_fasta()

# generates the PSSM
#gen.get_pssm(run=False)

# map the pssm to the pdb
#gen.map_pssm()

```

Run that on the cluster since the calculation of the PSSM files can take a few hours