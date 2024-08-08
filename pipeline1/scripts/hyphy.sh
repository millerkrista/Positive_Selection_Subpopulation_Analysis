#!/bin/sh

# Job Name
#$ -N hyphy-fubar

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

# run command: qsub -l avx2 scripts/hyphy.sh /hpcdata/vrc_vpds/millerkrs/pos_selection_subpopulation_structure_analysis/data/mexico fubar

# Set the directory containing the alignment and tree files
DATA_DIR=$1
# Set the hyphy method
METHOD=$2

module unload all
module load hyphy

# Iterate through each alignment file in the directory
for ALN in $DATA_DIR/*.aln; do
    # Get the corresponding tree file
    TRE="${ALN}.treefile"

    # Get the file name without the directory path
    FILENAME=$(basename "$ALN")

    # Run hyphy fubar on the alignment and tree file
    hyphy $METHOD --alignment "$ALN" --tree "$TRE" --code Universal 
done
