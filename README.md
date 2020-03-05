# PSSMGen

| fair-software.nl recommendations | Badges |
|:-|:-:|
| [1. Code Repository](https://fair-software.nl/recommendations/repository)       | [![GitHub URL](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/DeepRank/pssmgen) |
| &nbsp;                                                                          | [![GitHub](https://img.shields.io/github/last-commit/DeepRank/pssmgen)](https://github.com/DeepRank/pssmgen) |
| [2. License](https://fair-software.nl/recommendations/license)                  | [![License](https://img.shields.io/github/license/DeepRank/pssmgen)](https://github.com/citation-file-format/cff-converter-python) |
| [3. Community Registry](https://fair-software.nl/recommendations/registry)      | [![Research Software Directory](https://img.shields.io/badge/rsd-PSSMGen-00a3e3.svg&labelColor=gray&color=00a3e3)](https://research-software.nl/software/pssmgen) |
| &nbsp;                                                                          | [![PyPI](https://img.shields.io/pypi/v/pssmgen)](https://pypi.org/project/pssmgen/) |
| [4. Enable Citation](https://fair-software.nl/recommendations/citation)         | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3635711.svg)](https://doi.org/10.5281/zenodo.3635711) |
| [5. Code Quality Checklist](https://fair-software.nl/recommendations/checklist) | [![cii best practices](https://bestpractices.coreinfrastructure.org/projects/1811/badge)](https://bestpractices.coreinfrastructure.org/projects/1811)  |

<!--

| **Other**                                                                       | **Badge** |
| Continuous Integration                                                          | [![Build Status](https://travis-ci.org/research-software-directory/research-software-directory.svg?branch=master)](https://travis-ci.org/research-software-directory/research-software-directory) |
| &nbsp;                                                                          | [![Build status](https://ci.appveyor.com/api/projects/status/vki0xma8y7glpt09/branch/master?svg=true)](https://ci.appveyor.com/project/NLeSC/xenon-cli/branch/master)  |
| Code Analysis                                                                   | [![CodeClimate](https://api.codeclimate.com/v1/badges/ed3655f6056f89f5e107/maintainability)](https://codeclimate.com/github/DynaSlum/satsense/maintainability) |
| &nbsp;                                                                          | [![Codacy Badge](https://api.codacy.com/project/badge/Grade/6e3836750fe14f34ba85e26956e8ef10)](https://www.codacy.com/app/c-meijer/eEcoLiDAR?utm_source=www.github.com&amp;utm_medium=referral&amp;utm_content=eEcoLiDAR/eEcoLiDAR&amp;utm_campaign=Badge_Grade) |
| &nbsp;                                                                          | [![SonarCloud](https://sonarcloud.io/api/project_badges/measure?project=nlesc%3AXenon&metric=alert_status)](https://sonarcloud.io/dashboard?id=nlesc%3AXenon) |
| Code Coverage                                                                   | [![codecov](https://codecov.io/gh/wadpac/GGIR/branch/master/graph/badge.svg)](https://codecov.io/gh/wadpac/GGIR) |
| &nbsp; | [![SonarCloud](https://sonarcloud.io/api/project_badges/measure?project=xenon-middleware_xenon-grpc&metric=coverage)](https://sonarcloud.io/component_measures?id=xenon-middleware_xenon-grpc&metric=Coverage) |
| &nbsp; | [![Scrutinizer](https://scrutinizer-ci.com/g/NLeSC/mcfly/badges/coverage.png?b=master)](https://scrutinizer-ci.com/g/NLeSC/mcfly/statistics/) |
| &nbsp; | [![Coveralls](https://coveralls.io/repos/github/eEcoLiDAR/eEcoLiDAR/badge.svg)](https://coveralls.io/github/eEcoLiDAR/eEcoLiDAR) |
| &nbsp; | [![CodeClimate](https://api.codeclimate.com/v1/badges/ed3655f6056f89f5e107/test_coverage)](https://codeclimate.com/github/DynaSlum/satsense/test_coverage) |
| Documentation                                                                   | [![ReadTheDocs](https://readthedocs.org/projects/xenon-tutorial/badge/?version=latest)](https://xenon-tutorial.readthedocs.io/en/latest/?badge=latest) |

_(Customize these badges with your own links. Check https://shields.io/ to see which badges are available.)_

-->

Generates consistent PSSM and/or PDB files for protein-protein complexes

## Install

1. Make sure BLAST is installed and its database is available on your machine. Otherwise, install BLAST and download its databases by following the [BLAST guide](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
2. Install the PSSMgen by `pip install PSSMGen`.


## Requirements for file structures and names

`PSSMGen` is geared toward computing the pssm files for all models of a particular protein-protein complex.

### File structures
This tool assumes your files have following structure:

```
 workdir
 |_ pdb
 |_ fasta
 |_ pssm_raw
 |_ pssm
 |_ pdb_nonmatch
```

- `workdir` is your working directory for one specific protein-protein complex.
- `pdb` folder must contain the PDB files.
- `fasta` folder contains the protein sequence [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files. The code can generate the FASTA files by extracting sequences from the `pdb` file , or you can manually create this folder and put customed FASTA files there.
- `pssm_raw` folder stores the PSSM files. The code can automatically generate them, or you can manually create this folder and put customed PSSM files there.
- `pssm` folder stores consistent PSSM files, whose sequences are aligned with those of PDB files. This folder and its files are created automatically.
- `pdb_nonmatch` folder stores the inconsistent PDB files, while the related consistent PDB files are in the `pdb` folder. This folder and its files are created automatically.

### File names
The code assumes you follow the naming rules for different file types:
- PDB files:   caseID_*.chainID.pdb
- FASTA files: caseID.chainID.fasta
- PSSM files:  caseID.chainID.pssm, caseID_*.chainID.pdb.pssm


## Examples

### Get consistent PSSM and PDB files for a specific protein-protein complex

Here is an example for complex `7CEI`:

The file structure and input files should look like:
```
7CEI
└── pdb
    ├── 7CEI_1w.pdb
    ├── 7CEI_2w.pdb
    └── 7CEI_3w.pdb
```


```python
from pssmgen import PSSM

# initiate the PSSM object
gen = PSSM(work_dir='7CEI')

# set psiblast executable, database and other psiblast parameters (here shows the defaults)
gen.configure(blast_exe='/home/software/blast/bin/psiblast',
            database='/data/DBs/blast_dbs/nr_v20180204/nr',
            num_threads = 4, evalue=0.0001, comp_based_stats='T',
            max_target_seqs=2000, num_iterations=3, outfmt=7,
            save_each_pssm=True, save_pssm_after_last_round=True)

# generates FASTA files
gen.get_fasta(pdb_dir='pdb', chain=['A','B'], out_dir='fasta')

# generates PSSM
gen.get_pssm(fasta_dir='fasta', out_dir='pssm_raw', run=True)

# map PSSM and PDB to get consisitent files
gen.map_pssm(pssm_dir='pssm_raw', pdb_dir='pdb', out_dir='pssm', chain=['A','B'])

# write consistent files and move
gen.get_mapped_pdb(pdbpssm_dir='pssm', pdb_dir='pdb', pdbnonmatch_dir='pdb_nonmatch')
```

The code will automatically create `fasta`, `pssm_raw`, `pssm` and `pdb_nonmatch` folders and related files.

### Using customed FASTA files

You can use customed FASTA files intead of extracting them from PDB file.

The file structure and input files should look like
```
7CEI
└── fasta
    ├── 7CEI.A.fasta
    └── 7CEI.B.fasta
```


```python
from pssmgen import PSSM

# initiate the PSSM object
gen = PSSM('7CEI')

# set psiblast executable, database
gen.configure(blast_exe='/home/software/blast/bin/psiblast',
            database='/data/DBs/blast_dbs/nr_v20180204/nr')

# generates PSSM
gen.get_pssm()
```

### Using customed PSSM files

You can also use avaliable PSSM files intead of calculating them again.

The file structure and input files should look like
```
7CEI
├── pdb
│   ├── 7CEI_1w.pdb
│   ├── 7CEI_2w.pdb
│   └── 7CEI_3w.pdb
└── pssm_raw
    ├── 7CEI.A.pssm
    └── 7CEI.B.pssm
```

```python
from pssmgen import PSSM

# initiate the PSSM object
gen = PSSM('7CEI')

# map PSSM and PDB to get consisitent files
gen.map_pssm()

# write consistent files and move
gen.get_mapped_pdb()
```