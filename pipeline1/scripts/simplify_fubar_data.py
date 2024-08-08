#!/usr/bin/python3

"""Script that takes FUBAR data as input and outputs a text file containing descriptive headers and
all numeric data associated with codon site, dS, dN, and P(dS<dN)"""

import os
import json
import argparse


def main():
    # get command line arguments
    args = get_args()
   
    # gets json files from directory into list
    json_files = [f for f in os.listdir(args.directory) if f.endswith('.json')]
    
    # for each json file in list
    for file in json_files:
        
        # open file and load data
        json_data = load_json_data(f'{args.directory}/{file}')
        
        # write relevant data to new file
        write_data_to_file(file, args.directory, json_data)


def get_args():
    """Get command line arguments"""
    parser = argparse.ArgumentParser(description='Extract relevant data from FUBAR, output to text')
    
    parser.add_argument('-d', '--dir',
    type = str, 
    dest = 'directory', 
    help = 'Directory to read json files')
    
    args = parser.parse_args() 
    return args
    
    
def load_json_data(filename):
    """Loads data from json file"""
    
    # open json file
    json_file = open(filename)
    
    # load json data into object
    json_data = json.load(json_file)
    
    # close file
    json_file.close()

    return json_data


def write_data_to_file(filename, out_directory, json_data):
    """Simplifies json data and outputs to text file"""
    
    # get headers as object
    headers = ['site', 'dS', 'dN', 'P(dS<dN)', 'dN/dS']
    
    # get all numeric data as object
    data = json_data['MLE']['content']['0']
    
    # get output file name
    newfile = filename.replace('.json', '_simple.txt')

    # initialize a counter for codon position
    counter = 1
    
    # initialize tab character to use in f-string
    tab = '\t'
    
    # open output file
    with open(f'{out_directory}/{newfile}', 'a') as outfile:

        # write headers to file, tab-delimited
        outfile.write('\t'.join(str(header) for header in headers))
        outfile.write('\n')

        # for each row of data,
        for row in data:
            # get dN/dS
            dN_dS = row[1] / row[0]
            # get a subset that is dS, dN, and P(dS<dN)
            dat = [row[0], row[1], row[4], dN_dS]
            # write codon position and associated data to file, tab-delimited
            outfile.write(f"{counter}\t{tab.join(map(str, dat))}\n")
            # step counter
            counter += 1


if __name__ == "__main__":
    main()


