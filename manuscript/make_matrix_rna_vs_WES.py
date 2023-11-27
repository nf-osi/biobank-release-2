#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import gzip
import pandas as pd

__author__ = "Alexandra Scott (sasha.scott@sagebase.org)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2023-09-26 11:41 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
    	make_matrix_rna_vs_WES.py\n\
    	author: " + __author__ + "\n\
    	version: " + __version__ + "\n\
    	description: convert relatedness data from a line by line input to matrix format")
    parser.add_argument('-i', '--input',
                        metavar='FILE', dest='input_path',
                        type=str, default=None,
                        help='Input file containing tab separated sample IDs in first two columns and relatedness coefficient in the third column')

    # parse the arguments
    args = parser.parse_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)

    # send back the user input
    return args
    
# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data
    
    
def get_relatedness(relatedness_file):
    relatedness = {}
    for l in relatedness_file:
        l_split = l.rstrip('\n').split('\t')
        if l_split[0] not in relatedness:
            relatedness[l_split[0]] = {}
        relatedness[l_split[0]][l_split[1]]=l_split[2]
    r = pd.DataFrame(relatedness)
    return r

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        input_file = sys.stdin
    else:
        input_file = get_file(args.input_path)

    # run functions
    relatedness = get_relatedness(input_file) 
    print(relatedness.to_csv(na_rep="NA", sep='\t'))
 
    # close files
    input_file.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except OSError:
        if e.errno != 32:
        # ignore SIGPIPE
            raise
