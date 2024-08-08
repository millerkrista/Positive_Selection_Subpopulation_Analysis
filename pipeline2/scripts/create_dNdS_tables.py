#!/usr/bin/env python

"""Script to output files showing shared codon positions where positive selection is occurring
for a given set of participants at annotation sites with Tajima's D >= 1.5. There will be four output files
per annotation: *_dNdS.csv and *_probability.csv (comma separated text files for making heatmaps), 
*_dNdS_table.png, *_probability_table.png (image files for viewing tables)"""

import os
import re
import json
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt


def main():
    
    # get filenames as list for every consensus_annotations file
    files = [file for file in os.listdir('data/') if file.endswith('consensus_annotations.csv')]
    
    # use filenames to get individual IDs
    ids = [f.replace('_consensus_annotations.csv', '') for f in files]
    
    # from list of IDs, call the function to create a dataframe
    dNdS_df = create_selection_df(ids)
    prob_df = create_selection_df(ids)
    
    # for every ID, call function to fill in respective dataframe column
    for id in ids:
        dNdS_df = fill_selection_df(dNdS_df, id, dir='data', value='dN/dS')
        prob_df = fill_selection_df(prob_df, id, dir='data', value='P(dS<dN)')
    
    # save entire selection_df as csv
    dNdS_df.to_csv('results/tables/pos_selection_sites_dNdS.csv')
    prob_df.to_csv('results/tables/pos_selection_sites_probability.csv')
        
    # get dictionary of significant annotations and their positions
    sig_annotations = get_sig_annotations(dir='data')
    
    # go through ranges where significant annotation
    for annotation in sig_annotations.keys():
        
        # get subset of dNdS containing only range of annotation
        dNdS_table = dNdS_df[(dNdS_df.index >= sig_annotations[annotation][0]) & (dNdS_df.index <= sig_annotations[annotation][1])]
        # get subset of posterior probability containing only range of annotation
        prob_table = prob_df[(prob_df.index >= sig_annotations[annotation][0]) & (prob_df.index <= sig_annotations[annotation][1])]
        
        # replace / characters from annotation name before saving
        annotation_save = annotation.replace('/', '_')
        
        # save as csv
        dNdS_table.to_csv(f'results/tables/{annotation_save}_dNdS.csv')
        prob_table.to_csv(f'results/tables/{annotation_save}_probability.csv')
        
        if not dNdS_table.empty:
            # TODO: CHANGE WHERE IMAGES SAVED
            # save as table image
            dataframe_to_table_image(dNdS_table, annotation, suffix='_dNdS')
            dataframe_to_table_image(prob_table, annotation, suffix='_probability')

    
def create_selection_df(ids):
    """Function for creating an empty dataframe for a set of IDs"""
        
    # create an empty dataframe with participant IDs as columns
    selection_df = pd.DataFrame(columns = ids)
    
    return selection_df
    

def fill_selection_df(selection_df, id, dir, value):
    """Function to fill table with specified value for significant codons (their HXB2 reference codon equivalent) for a given participant"""
    
    # open the consensus annotations file associated with the ID
    consensus_annotations = pd.read_csv(f'{dir}/{id}_consensus_annotations.csv')
    
    # change index to start from 1
    consensus_annotations.index = consensus_annotations.index + 1
    
    # Create a regular expression pattern to match the fubar file name
    pattern = rf'{id}.*\.FUBAR_simple\.txt'

    # Get the list of files in the directory that match the pattern
    files = [file for file in os.listdir(dir) if re.match(pattern, file)]

    # Read the first matching file
    fubar_data = pd.read_table(os.path.join(dir, files[0]))

    # Get the codon sites as a list
    sites = fubar_data['site'].tolist()

    for site in sites:
        # Get index where aligned consensus matches site
        site_index, og_consensus_decimal = get_site_index(site, consensus_annotations)

        # Get the associated hxb2 codon position
        hxb2_pos = get_hxb2_position(site_index, og_consensus_decimal, consensus_annotations)

        # Add value to correct hxb2 position and participant ID of selection dataframe
        val = fubar_data.loc[fubar_data['site'] == site, value].values[0]
        selection_df.at[hxb2_pos, id] = val

    return selection_df.sort_index()


def get_site_index(site, consensus_annotations):
    """Function takes as input a sequence site, checks what codon is at that site in the original consensus,
    and outputs the index where that codon is located in the aligned consensus. If the position in the original consensus is a gap,
    function will output the index of the left/lower 'real' codon and a decimal to represent what position the gap is between the left and right
    'real' consensus codons."""
    
    # see what original consensus codon is located at the site
    og_consensus_pos = consensus_annotations.at[site, 'Original.Consensus.Position']
    
    
    # original consensus position is a gap
    if np.isnan(og_consensus_pos):
        
        # find the nearest lower, upper index where Original Consensus contains a 'real' codon
        nearest = (consensus_annotations.loc[:site, 'Original.Consensus.Position'].last_valid_index(), 
                   consensus_annotations.loc[site:, 'Original.Consensus.Position'].first_valid_index())
        
        lower = nearest[0]
        upper = nearest[1]
        
        
        if lower is None:
            # significant site in gap at beginning
            # set decimal to -0.5
            og_consensus_decimal =  -0.5
            
        
        elif upper is None:
            # significant site in gap at end
            # set decimal to 0.5
            og_consensus_decimal =  0.5
            
                    
        else:
            # this is the decimal between the two 'real' codons that our significant site is located at
            og_consensus_decimal = ((upper - lower) - (upper - site)) / (upper - lower)
    
        
        # take the left/lower real codon as the original consensus position
        og_consensus_pos = consensus_annotations.at[lower, 'Original.Consensus.Position']
    
        # get the index where the original consensus codon exists in the aligned consensus
        site_index = consensus_annotations.index[consensus_annotations['Consensus.Aligned.to.HXB2.Position'] == og_consensus_pos][0]
    
    
    else:
        # original consensus position is not a gap
        # set decimal to 0
        og_consensus_decimal = 0
        
        # get the index position where the codon exists in the aligned consensus
        site_index = consensus_annotations.index[consensus_annotations['Consensus.Aligned.to.HXB2.Position'] == og_consensus_pos][0]
    
    return site_index, og_consensus_decimal
    
    
def get_hxb2_position(site_index, og_consensus_decimal, consensus_annotations):
    """Function to get the HXB2 codon position where significant dN/dS was found;
    if position is a gap, make position a decimal between the nearest 'real' codons""" 
    
    # get hxb2 position from the site index
    hxb2_pos = consensus_annotations.at[site_index, 'HXB2.Position']
    
    # hxb2_pos is a gap
    if np.isnan(hxb2_pos):
        
        # find the nearest upper, lower index where HXB2.Position contains a 'real' codon
        nearest = (consensus_annotations.loc[:site_index, 'HXB2.Position'].last_valid_index(), 
                   consensus_annotations.loc[site_index:, 'HXB2.Position'].first_valid_index())
                   
        lower = nearest[0]
        upper = nearest[1]
        
        
        if lower is None:
            # significant codon is before hxb2 start
            # subtract 0.5 from the right/upper 'real' codon to get hxb2 position
            hxb2_pos = consensus_annotations.at[upper, 'HXB2.Position'] - 0.5
            
        elif upper is None:
            # significant codon occurs after hxb2 end
            # add 0.5 to the left/lower 'real' hxb2 codon to get hxb2 position
            hxb2_pos = consensus_annotations.at[lower, 'HXB2.Position'] + 0.5
                    
        else:
            # this is the decimal between the two 'real' codons that our significant site is located at
            hxb2_decimal = ((upper - lower) - (upper - site_index)) / (upper - lower)
            
            # add the decimal to the left/lower 'real' hxb2 codon to get hxb2 position
            hxb2_pos = consensus_annotations.at[lower, 'HXB2.Position'] + hxb2_decimal
            
        
        # return hxb2_pos
        return hxb2_pos
        
    
    else:
        # return hxb2 position if it is not a gap
        # add consensus fraction to the position
        return hxb2_pos + og_consensus_decimal
        
        
def get_sig_annotations(dir):
    """Function to get dictionary of annotation and its HXB2 start/stop codon position when the region had a significant Tajima's D;
    input is the selection dataframe, returns a dictionary with format {annotation: [start, stop]}"""
    
    with open('results/sig_midpoints.json', 'r') as json_file:
        # load json file into dictionary
        sig_midpoints_dict = json.load(json_file)
    
    
    # initialize dictionary to hold annotation and positions
    sig_annotations_dict = {}
    

    for id in sig_midpoints_dict.keys():
        
        # open the consensus annotations file for the ID
        consensus_annotations = pd.read_csv(f'{dir}/{id}_consensus_annotations.csv')
        
        for midpoint in sig_midpoints_dict[id]:
            
            # get codon position of midpoint by dividing by 3, rounding down to nearest whole number
            midpoint = midpoint // 3
            
            # get midpoint +- 17 to get the whole Tajimas D window
            D_start = max(1, midpoint - 17)  # if less than 0, set to 1
            D_stop = midpoint + 17
            
            # get indices where Tajimas D window is located in the aligned consensus
            D_start_idx, start_decimal = get_site_index(D_start, consensus_annotations)
            D_stop_idx, stop_decimal = get_site_index(D_stop, consensus_annotations)
            D_start_idx = math.floor(D_start_idx + start_decimal)
            D_stop_idx = math.ceil(D_stop_idx + stop_decimal)
            
            # subset dataframe to just show annotation regions associated with significant Tajimas D window
            ann_df = consensus_annotations.loc[D_start_idx: D_stop_idx, ['Regions1', 'Regions2', 'Regions3']]
            
            # initialize list to hold annotations associated with Tajimas D range
            annotations = []
            
            # get real annotations in Regions1, 2, and 3 as a list of lists
            for region in ann_df:
                annotations.append(ann_df[region].unique().tolist())
            
            # flatten list, remove NAs
            annotations = [a for lis in annotations for a in lis if not pd.isnull(a)]
            
            
            # for each annotation
            for annotation in annotations:
                
                if annotation in sig_annotations_dict.keys():
                    # if annotation already in dict, move to next annotation
                    continue
                
                else:
                    # get positions of annotation as values
                    ann_pos = consensus_annotations.loc[consensus_annotations['Regions1'].isin([annotation]) | consensus_annotations['Regions2'].isin([annotation]) | consensus_annotations['Regions3'].isin([annotation]), 'HXB2.Position'].tolist()
                    
                    # add to dictionary with format {annotation: [start position, stop position]}
                    sig_annotations_dict[annotation] = [int(ann_pos[0]), int(ann_pos[-1])]
        
            
            # check significant annotations to handle edge cases (beginning/end of hxb2)
            for annotation in sig_annotations_dict.keys():
                
                if sig_annotations_dict[annotation][0] == 1:
                    # change annotation region to start from 0 if it is beginning of sequence (i.e signal peptide will start at 0)
                    sig_annotations_dict[annotation][0] = 0
        
                if sig_annotations_dict[annotation][1] == consensus_annotations.at[consensus_annotations.loc[:, 'HXB2.Position'].last_valid_index(), 'HXB2.Position']:
                    # change annotation region to end at the end of consensus sequence if it is the end of hxb2 sequence
                    sig_annotations_dict[annotation][1] = len(consensus_annotations['Original.Consensus.Position'])
                    
    
    # return dictionary of annotations and ranges
    return sig_annotations_dict   


def dataframe_to_table_image(table_df, annotation, suffix):
    """Function takes as input a dictonary of annotations and their start/stop codon positions;
    creates a subset of the selection dataframe and saves subset as a table image"""
    
    # create a table from dataframe
    fig, ax = plt.subplots(figsize=(45, 20))  # Adjust the figure size as needed
    ax.axis('off')  # Remove axis
    ax.set_title(annotation)
   
    # replace NaN cells in table with whitespace
    table_df = table_df.replace(np.NaN, '')
    
    table = pd.plotting.table(ax, table_df, loc='center', cellLoc='center')

    # Adjust the font size of the table
    table.auto_set_font_size(False)
    table.set_fontsize(8)  # Adjust the font size as needed

    # Adjust the cell size
    table.scale(1.5, 2)  # Increase the cell size by a factor of 1.5 (adjust as needed)

    # replace / characters from annotation name before saving
    annotation = annotation.replace('/', '_')

    plt.show()
    plt.savefig(f'results/plots/{annotation}{suffix}_table.png', bbox_inches = "tight")
    plt.close()
    
    
if __name__ == "__main__":
    main()

