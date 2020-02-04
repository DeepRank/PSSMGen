import os
from setuptools import setup

cwd = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit pssmgen/__version__.py
version = {}
with open(os.path.join(cwd, 'pssmgen', '__version__.py')) as f:
    exec(f.read(), version)

setup(
    name='PSSMGen',
    description='Generate PSSM and matched PDB files for protein-protein complexes'
    version=version['__version__'],
    url='https://github.com/DeepRank/PSSMGen',
    packages=['pssmgen'],

    install_requires=[
        'numpy >= 1.13',
        'scipy',
        'biopython',
        'pdb2sql'],

    extras_require= {
        'test': ['nose', 'coverage', 'pytest',
                 'pytest-cov','codacy-coverage','coveralls'],
    }
)
