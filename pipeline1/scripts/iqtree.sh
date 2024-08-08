#!/bin/sh

# Job Name
#$ -N iqtree

# Execute the script from the current working directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

#$ -S /bin/bash

# Send the output of the script to a directory called 'UGE-output' in the current working directory (cwd)
#$ -o UGE-output/

# Send mail when the job is submitted, and when the job completes
# -m be

#  Specify an email address to use
#$ -M millerkrs@nih.gov

# run command: qsub scripts/iqtree.sh /hpcdata/vrc_vpds/millerkrs/pos_selection_subpopulation_structure_analysis/data/VRC601

# Set the directory containing the alignment and tree files
DATA_DIR=$1

module unload all
module load iqtree

# Iterate through each alignment file in the directory
for ALN in $DATA_DIR/*.aln; do

    # Run iqtree on the alignment
    iqtree -s $ALN -m MFP –bb 1000 -nt AUTO 

done
