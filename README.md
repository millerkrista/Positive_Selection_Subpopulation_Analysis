# Positive Selection and Subpopulation Analysis of HIV-1 *env* gene
This project was done as part of the National Institutes of Health's Graduate Data Science Summer Program. I worked under Frida Belinky, a computational biologist with Eli Boritz's lab at the National Institute of Allergy and Infectious Disease. With Frida's mentorship, I created a novel workflow for analyzing high-throughput single-genome sequencing data from 40 participants with chronic HIV infection. 

## Background
When looking at sequence alignments of a single participant's HIV-1 population, there seemed to be distinct subpopulations formed by the presence of different haplotypes in *env* (image of example alignment shown below).

<img width="395" alt="Screenshot 2025-02-21 at 12 57 04 PM" src="https://github.com/user-attachments/assets/72507ed2-6a9a-44cd-991a-37464e22981d" />

We had some ideas about what might be driving host populations to evolve into these subpopulations: 

1. Positive Selection: Mutations to DNA that alter the protein are maintained in a population, perhaps because they increase fitness in some way. In our case, a fitness increase could refer to improved capability of avoiding host antibodies or improved binding to coreceptors important for infection.

2. Balancing Selection: Multiple haplotypes are maintained in a population, perhaps because they are all adaptive in some way. In our case, one combination of amino acids that make up a haplotype could improve the virus's ability to escape host antibodies while another increases coreceptor binding affinity.
