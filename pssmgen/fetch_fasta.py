#!/usr/bin/env python

"""
Get the official FASTA file from PDB website for each chain of the given PDB ID.

Usage: python fetch_fasta.py <PDB IDs>
Example: python fetch_fasta.py 1CBH 7CEI 4CPA

Author: {0} ({1})
"""

__author__ = "Cunliang Geng"
__email__ = "gengcunliang AT gmail.com"
USAGE = __doc__.format(__author__, __email__)

import sys
import re
import requests


def check_input(args):
    """Validate user input
    
    Arguments:
        args {tuple} --  user input arguments
    """
    if not len(args):
        sys.exit(USAGE)


def write_chain_sequence(pdbid):
    """Write the complete sequence of each chain to file for a PDB ID.

    Args:
        pdbid: Official ID of PDB file in Protein Data Bank.
    Output:
        pdbID.chainID.fasta files.
    """

    # check the format of PDB ID
    rule = re.compile("[0-9a-zA-Z]{4}")

    if rule.match(pdbid):
        url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList={}&compressionType=uncompressed".format(pdbid)
        r = requests.get(url)

        # check if the returned content is fasta or not
        if r.content.startswith(b">"):
            for i in r.content.decode().split(">")[1:]:
                f = i[:4] + "." + i[5] + ".fasta"
                fh = open(f, "w")
                fh.write(">" + i)
                fh.close()
                print(f)

        else:
            try:
                raise NameError("Warning: The PDB ID '{}' currently NOT exist in the Protein Data Bank.".format(pdbid))
            except NameError as e:
                print(e)
    else:
        try:
            raise NameError("Warning: Wrong format of PDB ID '{}'. It should be a 4-character identifier with only alphabet and/or number.".format(pdbid))
        except NameError as e:
            print(e)


def write_chain_sequence_multiPDB(*pdbids):
    """Write the complete sequence of each chain to file for multiple PDB IDs.

    Args:
        pdbids: a list of pdbIDs.
    """

    for i in pdbids:
        write_chain_sequence(i)


if __name__ == "__main__":
    check_input(sys.argv[1:])
    # print(sys.argv[1:])
    write_chain_sequence_multiPDB(*sys.argv[1:])