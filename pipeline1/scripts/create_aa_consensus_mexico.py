#!/usr/bin/env python3

"""Program to create consensus sequences for every protein alignment file in directory, outputs a file containing all consensus sequences for experiment"""

import os
import argparse
import numpy as np
import pandas as pd
import argparse
import shutil



def main():
    # get command line arguments
    args = get_args()
    
    # load in alignment files as list
    aln_files = [f for f in os.listdir(args.directory) if f.endswith('.aln.translated')]

    
    # for each file in the list
    for file in aln_files:
        
        individual_id = file.replace('.aln.translated', '')
        
        # get list of sequences
        sequences = get_sequences_from_file(f'{args.directory}/{file}')
        
        # get amino acid frequency counts
        counts = get_counts_matrix(sequences)
        
        # reconstruct sequence from count matrix
        most_common = counts[counts == counts.max(axis=0)]
        consensus = most_common.apply(counts_to_sequence, axis=0)

        # make copy of alignment file
        shutil.copyfile(f'{args.directory}/{individual_id}.aln.translated', f'{args.directory}/{individual_id}.aln.translated.consensus')
        
        # write consensus sequence to copied file
        with open(f'{args.directory}/{individual_id}.aln.translated.consensus', 'a') as outfile1:
            outfile1.write(f'>{individual_id} Consensus Sequence\n')
            outfile1.write(''.join(consensus))
        with open(f'{args.directory}/{args.outfile}', 'a') as outfile2:
            outfile2.write(f'\n>{individual_id} Consensus Sequence\n')
            outfile2.write(''.join(consensus))


def get_args():
    """Get command line arguments"""
    parser = argparse.ArgumentParser(description='Create consensus sequences from MSA files based on amino acid frequency')
    
    parser.add_argument('-o', '--outfile', dest='outfile', help='File to store consensus sequences')
    
    parser.add_argument('-d', '--directory', dest='directory', help='Directory to read alignment files')
    
    args = parser.parse_args()
    
    return args
    
    
def get_counts_matrix(sequences):
    """From a list of DNA sequences, creates a matrix counting amino acid frequencies at each position of a DNA sequence;
    @input: list of sequences from MSA
    @return: amino acid frequency matrix at each position"""
    
    # initialize a matrix of zeros with length = sequence length
    count_matrix = np.zeros([22, len(sequences[0])], dtype = int)
    
    # create a dataframe from zeros matrix with amino acid options as indices
    counts = pd.DataFrame(count_matrix, index = ['a', 'c', 'd', 'e', 'f', 'g', 'h', 
    'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y', 'x', '-'])
    
    # for each sequence from MSA
    for seq in sequences:
        
        # look at each amino acid in sequence
        for ind, aa in enumerate(seq):
            
            # add 1 to counts matrix for observed amino acid
            counts[ind][aa] += 1
            
    return counts
    
    
def get_sequences_from_file(aln_file):
    """From a MSA file, creates a list of sequences
    @input: name of alignment file
    @return: list of DNA sequences"""
    
    # open file
    with open(aln_file, 'r') as aln_file:
        
        # initialize list of sequences for file
        sequences = []
        
        # read each line of file
        for line in aln_file.readlines():
            
            # skip header lines
            if line.startswith('>'):
                continue
            
            # append sequences to sequence list
            else: 
                sequences.append(line.lower().strip())
        
        return sequences


def counts_to_sequence(pos):
    top_aa = pos.dropna().index
    key = ''.join(sorted(top_aa))
    return key


if __name__ == "__main__":
    main()


