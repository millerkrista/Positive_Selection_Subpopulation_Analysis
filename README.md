# Positive Selection and Subpopulation Analysis of HIV-1 *env* gene
This project was done as part of the National Institutes of Health's Graduate Data Science Summer Program. I worked under Frida Belinky, a computational biologist with Eli Boritz's lab at the National Institute of Allergy and Infectious Disease. With Frida's mentorship, I created a novel workflow for analyzing high-throughput single-genome sequencing data from 40 participants with chronic HIV infection. 

## Background
When looking at sequence alignments of a single participant's HIV-1 population, there seemed to be distinct subpopulations formed by the presence of different haplotypes in *env* (image of example alignment shown below).

<p align="center">
  <img width="395" alt="Screenshot 2025-02-21 at 12 57 04 PM" src="https://github.com/user-attachments/assets/72507ed2-6a9a-44cd-991a-37464e22981d" />
</p>

We had some ideas about what might be driving host populations to evolve into these subpopulations: 

1. Positive Selection: Mutations to DNA that alter the protein are maintained in a population, perhaps because they increase fitness in some way. In our case, a fitness increase could refer to improved capability of avoiding host antibodies or improved binding to coreceptors important for infection.

2. Balancing Selection: Multiple haplotypes are maintained in a population, perhaps because they are all adaptive in some way. In our case, one combination of amino acids that make up a haplotype could improve the virus's ability to escape host antibodies while another increases coreceptor binding affinity.

## Biological Question
We wanted to investigate these subpopulations to see where they were occurring and what could be causing them. 

## Process
The general steps in image below were performed by others, I started out with sequence alignments for each participant containing just their *env* sequences.

<p align="center">
  <img width="362" alt="Screenshot 2025-02-21 at 1 25 19 PM" src="https://github.com/user-attachments/assets/1c2c2707-b9e7-4463-98fd-82ec8fa79940" />
</p>

### Preliminary
Used DNASP to calculate the Tajima's D values on a sliding window of 100bp for each participant alignment file. Tajima's D is typically used to measure whether balancing selection is occurring. 

Used HYPHY FUBAR to calculate the dN/dS, or the ratio of nonsynonymous mutation to synonymous mutation rate. This is typically used to infer whether positive selection is occurring at a given site. 

### Pipeline 1
Tools/Languages/Dependencies: Python, R

This workflow takes as input the positive selection data we got from HYPHY FUBAR and the participant alignments. We translate the alignments to protein and obtain consensus sequences for each participant, align the consensus sequences to the reference HXB2, and output files containing codon position mappings for each participant alignment compared to reference HXB2, along with relevant gene annotations for those positions. These files will be utilized later.

This workflow also creates plots for each participant showing dN/dS and Tajima's D values across *env* (example image shown below). These plots were intended to indicate areas where positive selection might be driving subpopulation formation.

<p align="center">
  <img width="653" alt="Screenshot 2025-02-21 at 2 38 47 PM" src="https://github.com/user-attachments/assets/e04d6105-7a29-4d78-93ab-90de979d5d55" />
</p>


### Pipeline 2
Tools/Languages/Dependencies: Python, R

This workflow takes as input the Tajima's D data we obtained in the preliminary stage and the position/annotation mapping files we created in pipeline1. 

### Pipeline 3



