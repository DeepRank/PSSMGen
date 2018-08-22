#from distutils.core import setup
from setuptools import setup

setup(
    name='PSSMGen',
    description='Generates PSSM files for deep learning applications of proteibn-protein docking',
    version='0.0',
    url='https://github.com/DeepRank/PSSMGen',
    packages=['pssmgen'],


    install_requires=[
        'numpy >= 1.13',
        'scipy'],
        #'tarfiles',
        #'pickle'],

    extras_require= {
        'test': ['nose', 'coverage', 'pytest',
                 'pytest-cov','codacy-coverage','coveralls'],
    }
)
