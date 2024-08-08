#!/usr/bin/env python
"""Script to translate nucleotide sequences in fasta file to protein and write to outfile"""

import os
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import Gapped



def main():
    
    # get command line arguments
    args = get_args()
    
    # read fasta file names into list
    aln_files = [f for f in os.listdir(args.directory) if f.endswith('_rev2miss.aln')]
    
    for file in aln_files:
        
        # read in sequences
        seqs = get_sequence_dict_from_file(f'{args.directory}/{file}')
        
        # get individual ID from file name
        individual_id = file.replace('_rev2miss.aln', '')
        
        outfile = f'{args.directory}/{individual_id}_rev2miss.aln.translated'
        
        # write translated sequences to outfile
        write_translated_seqs_to_file(outfile, seqs)


def get_args():
    """Get command line arguments"""
    parser = argparse.ArgumentParser(description='Translate fasta files to protein sequence')
    
    parser.add_argument('-d', '--dir',
    type = str, 
    dest = 'directory', 
    help = 'Directory to read alignment files')
    
    args = parser.parse_args() 
    return args
    

def get_sequence_dict_from_file(filename):
    '''Function to open fasta file and create dictionary with format {header: sequence....}'''
    
    # open file
    with open(filename, 'r') as seqfile:
        
        # get dictionary
        seq_dict = SeqIO.to_dict(SeqIO.parse(seqfile, "fasta"))
        
        return seq_dict
        

def write_translated_seqs_to_file(outfile, sequence_dict):
    
    # open file
    with open(outfile, 'a') as outfile:
        
        for header, seq in sequence_dict.items():
            # translate sequence
            seq = seq.translate(gap = '-')
            # write header and translated sequence to file
            outfile.write(f'>{header}\n{seq.seq}\n')
            

if __name__ == "__main__":
    main()

