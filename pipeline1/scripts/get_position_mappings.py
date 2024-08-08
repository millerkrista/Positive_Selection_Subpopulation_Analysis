#!/usr/bin/env python3

'''Script takes as input a file containing consensus sequences for each participant from experimental and 
a file containing those consensus sequences aligned to HXB2, returns files for each individual 
showing codon position mappings for HXB2, the aligned consensus, and the original consensus

Run with 'python3 get_seq_pos_mappings_from_alignments.py -c [file with consensus seqs] -a [file of consensus seqs aligned to HXB2]' '''

import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO


def main():
    # get command-line arguments
    args = get_args()
    
    # open file of consensus sequences created from alignments, get dictionary
    og_cons_dict = get_sequence_dict_from_file(args.consensus_file)
    
    # open file of consensus sequences aligned with hxb2, get dictionary
    cons_aln_with_hxb2_dict = get_sequence_dict_from_file(args.aln_file)
    
    # remove hxb2 from dict, get aligned hxb2 sequence as object
    hxb2_aln_seq = cons_aln_with_hxb2_dict.pop('MH758564.1')
    
    # for each individual identifier
    for id in og_cons_dict.keys():
        
        # get sequences associated with id
        aln_con_seq = cons_aln_with_hxb2_dict[id]
        og_con_seq = og_cons_dict[id]
        
        # get pandas dataframe that has position mappings
        pos_mapping_df = create_pos_mapping_df(hxb2_aln_seq, aln_con_seq, og_con_seq)
        
        # save that mapping to a file for later use in R
        pos_mapping_df.to_csv(f'{args.directory}/{id}_position_mappings.csv')
        
   
def get_args():
    """Get command line arguments"""
    parser = argparse.ArgumentParser(description='Map codon positions for participants to reference')
    
    parser.add_argument('-c', '--confile', 
    dest = 'consensus_file', 
    help = 'File to read consensus sequences, should be protein sequence')
    
    parser.add_argument('-a', '--alnfile', 
    dest = 'aln_file', 
    help = 'File to read HXB2 aligned consensus sequences, should be protein sequence')
    
    parser.add_argument('-d', '--dir', type = str, dest = 'directory', help = 'Directory to write position mapping files')
    
    args = parser.parse_args() 
    return args
    
    
def get_sequence_dict_from_file(filename):
    '''Function to open fasta file and create dictionary with format {header: sequence....}'''
   
    # open file
    with open(filename, 'r') as seqfile:
        
        # get dictionary
        seq_dict = SeqIO.to_dict(SeqIO.parse(seqfile, "fasta"))
        
        return seq_dict
        

def get_seq_positions(seq):
    '''Function to iterate through a protein sequence (that may contain gaps) and create an array
    that maps index position (of the array) to protein position of the sequence if it did not have gaps'''
    
    # remove newline characters from sequence, uppercase sequence
    seq = seq.upper()
    
    # initialize amino acid list
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y','X', 'B', 'Z', 'J']
    
    # initialize position counter
    pos_count = 0
    
    # create an array of length seq, all zeros
    seq_pos = np.zeros(len(seq))
    
    # for each position in sequence,
    for ind, aa in enumerate(seq):
       
        # if position is an amino acid,
        if aa in amino_acids:
            
            # add to count,
            pos_count += 1 
            
            # put new count into array
            seq_pos[ind] = pos_count
        
        else:
            # position is a gap ('-'), put NaN into array
            seq_pos[ind] = 'NaN'
    
    # returns array of sequence positions
    return seq_pos

    
def create_pos_mapping_df(hxb2_seq, aln_con_seq, og_con_seq):
    '''Function that takes as input the aligned hxb2, consensus seq aligned to hxb2, 
    and the consensus seq created from the alignment; returns a dataframe mapping codon positions of
    all sequences to each other'''
    
    # sends each seq to get array of amino acid positions
    hxb2_pos = get_seq_positions(hxb2_seq)
    aln_con_pos = get_seq_positions(aln_con_seq)
    og_con_pos = get_seq_positions(og_con_seq)
    
    # check which sequence is longest
    l = [hxb2_pos, aln_con_pos, og_con_pos]
    length = len(max(l, key = len))
   
    # creating a list of index and column names 
    index_values = ['HXB2 Position', 'Consensus Aligned to HXB2 Position', 'Original Consensus Position'] 
    # use the length value created earlier
    column_values = range(1, length + 1)  
  
    # creating the dataframe 
    pos_mapping_df = pd.DataFrame(data = [hxb2_pos, aln_con_pos, og_con_pos], columns = column_values, index = index_values) 
    # transpose the dataframe, columns become rows and rows become columns
    pos_mapping_df = pos_mapping_df.transpose()
    
    # returns finished dataframe
    return pos_mapping_df
    

if __name__ == "__main__":
    main()


