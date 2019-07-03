#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Split a file into subsets based on the number of lines
"""
import sys
args = sys.argv[1:]
inputFile = args[0] # "targets.bed"

# number of output files desired
num_chunks = int(args[1]) # 5

# need to divide the file this many times
chunk_divisor = num_chunks - 1

# doesnt work with less than 1 chunk
if num_chunks < 1:
    raise

def num_lines(fname):
    """
    Gets the total number of lines in the file
    """
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return( i + 1 )

# total number of file lines
inputFile_numlines = num_lines(inputFile)

# the size that each chunk should be, except for the remainder
chunk_size  = inputFile_numlines // chunk_divisor

# read the input lines
with open(inputFile) as f:
    for i, line in enumerate(f):
        # chunk number = current line num / even size of chunks + 1
        chunk_num = ( i // chunk_size ) + 1

        # write line to output file
        output_file = "{0}.{1}".format(inputFile, chunk_num)
        with open(output_file, "a") as fout:
            fout.write(line)
