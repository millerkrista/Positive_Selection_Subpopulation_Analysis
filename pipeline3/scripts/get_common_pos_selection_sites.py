#!/usr/bin/env python

"""Script to go through table of P(dS<dN) values for entire gene, extract rows (codon positions) where 
given number of participants or more had significant dN/dS"""


# imports
import argparse
import pandas as pd


def main():
    # get CLI arguments
    args = get_args()
    
    # open file as dataframe
    df = pd.read_csv(args.infile, index_col = 0)
        
    # create a boolean mask where values are True if they are greater than or equal to 0.9
    mask = df >= 0.9
        
    # Sum the boolean values across columns (patients) for each row
    row_sums = mask.sum(axis=1)
        
    # define the minimum number of patients that should exceed 0.9
    min_patients = args.threshold
        
    # filter dataframe by rows where [threshold] columns or more have P(dS<dN) > 0.9
    filtered_df = df[row_sums >= min_patients]
    
    # save filtered dataframe as table
    filtered_df.to_csv(args.outfile)
    

def get_args():
    """parse CLI arguments"""
    
    parser = argparse.ArgumentParser(description='Create table of common significant positive selection sites')
    
    parser.add_argument('-i', '--infile', dest='infile', help='File containing P(dS<dN) values for gene')
    parser.add_argument('-t', '--threshold', dest='threshold', help='Threshold for how many participants need significance for codon to be extracted', type=int, default=5)
    parser.add_argument('-o', '--outfile', dest='outfile', help='File to write results')
    
    args = parser.parse_args()
    
    return args
    

if __name__ == "__main__":
    main()
    
