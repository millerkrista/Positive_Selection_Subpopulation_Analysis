# Positive Selection and Subpopulation Analysis of HIV-1 *env* Gene in Participants with Chronic Infection
This project was done as part of the National Institutes of Health Intramural Program's Graduate Data Science Summer Program. I worked under Frida Belinky, a computational biologist with Eli Boritz's lab at the National Institute of Allergy and Infectious Disease. With Frida's invaluable mentorship, I created a novel workflow for analyzing high-throughput single-genome sequencing data from 40 participants with chronic HIV infection. 


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

### Preliminary Steps
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

This workflow takes as input the Tajima's D data we obtained in the preliminary stage and the position/annotation mapping files we created in pipeline1. The output will be tables and heatmaps showing the probability of positive selection (P(dS<dN)) in areas of the sequence where at least one participant had a Tajima's D value greater than 1.5 (first example heatmap below). I also created a heatmap showing the probability of positive selection throughout all of *env* for every participant (second example heatmap below). 

These heatmaps were used to show that significant positive selection is occurring throughout *env* and that many of these regions with positive selection also have subpopulations forming. 

<p align="center">
  <img width="794" alt="Screenshot 2025-02-21 at 2 48 57 PM" src="https://github.com/user-attachments/assets/56a21d3f-c8a1-4bab-9ed3-32e11ebf904c" />
</p>

<p align="center">
  <img width="153" alt="Screenshot 2025-02-21 at 2 51 15 PM" src="https://github.com/user-attachments/assets/252b4fb9-c859-4d2b-b153-b68dde19c834" />
</p>

### Pipeline 3
Tools/Languages/Dependencies: Python

This workflow takes as input the file dN/dS values throughout *env*, a threshold of participants that should have positive selection occurring, and an outfile destination. The output will be a file showing codon positions where a threshold number of participants or more had significant positive selection occurring. 

These results were utilized in literature review--we were investigating whether these particular *env* positions had previously been annotated as useful for antibody escape, association with coreceptor binding, or any other functional property that would indicate why that position would be under significant positive selection in multiple participants with chronic HIV infection.


## Conclusions
More investigation required than what could be done in a three month summer project. It seemed that some subpopulations found in some participants were associated with coreceptor usage (CCR5 vs CXCR4), which was investigated using WEBPSSM and phylogenetic tree analysis. Other positions seemed to be associated with antibody escape. 

Research is ongoing at the Vaccine Research Center.


## Acknowledgements
<img width="176" alt="Screenshot 2025-02-21 at 3 07 25 PM" src="https://github.com/user-attachments/assets/7c5fb3d4-1a80-4376-88d8-89b218975586" />

<img width="122" alt="Screenshot 2025-02-21 at 3 06 51 PM" src="https://github.com/user-attachments/assets/b881dbfe-04b1-4c4f-9298-635ada884c26" />

The Virus Persistence and Dynamics Section at the Vaccine Research Center
- Eli Boritz
- Frida Belinky
- Pierce Radecki
- Sung Hee Ko
- Vanessa Guerra Canedo
- Ayushman Dobhal
- Kaelo Seatla
- Prakriti Mudvari
- Liliana Perez Rodriguez
- Pavitra Ramdas
- Marc Theberge


### References
Arrildt KT, Joseph SB, Swanstrom R. The HIV-1 env protein: a coat of many colors. Curr HIV/AIDS Rep. 2012 Mar;9(1):52-63. doi: 10.1007/s11904-011-0107-3. PMID: 22237899; PMCID: PMC3658113.

Yandrapally S, Mohareer K, Arekuti G, Vadankula GR, Banerjee S. HIV co-receptor-tropism: cellular and molecular events behind the enigmatic co-receptor switching. Crit Rev Microbiol. 2021 Aug;47(4):499-516. doi: 10.1080/1040841X.2021.1902941. Epub 2021 Apr 26. PMID: 33900141.

Ko SH, Bayat Mokhtari E, Mudvari P, Stein S, Stringham CD, Wagner D, Ramelli S, Ramos-Benitez MJ, Strich JR, Davey RT Jr, Zhou T, Misasi J, Kwong PD, Chertow DS, Sullivan NJ, Boritz EA. High-throughput, single-copy sequencing reveals SARS-CoV-2 spike variants coincident with mounting humoral immunity during acute COVID-19. PLoS Pathog. 2021 Apr 8;17(4):e1009431. doi: 10.1371/journal.ppat.1009431. PMID: 33831133; PMCID: PMC8031304.

Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013 Apr;30(4):772-80. doi: 10.1093/molbev/mst010. Epub 2013 Jan 16. PMID: 23329690; PMCID: PMC3603318.

Murrell B, Moola S, Mabona A, Weighill T, Sheward D, Kosakovsky Pond SL, Scheffler K. FUBAR: a fast, unconstrained bayesian approximation for inferring selection. Mol Biol Evol. 2013 May;30(5):1196-205. doi: 10.1093/molbev/mst030. Epub 2013 Feb 18. PMID: 23420840; PMCID: PMC3670733.

Rozas J, Ferrer-Mata A, Sánchez-DelBarrio JC, Guirao-Rico S, Librado P, Ramos-Onsins SE, Sánchez-Gracia A. DnaSP 6: DNA Sequence Polymorphism Analysis of Large Data Sets. Mol Biol Evol. 2017 Dec 1;34(12):3299-3302. doi: 10.1093/molbev/msx248. PMID: 29029172.

Jensen MA, Li FS, van 't Wout AB, Nickle DC, Shriner D, He HX, McLaughlin S, Shankarappa R, Margolick JB, Mullins JI. Improved coreceptor usage prediction and genotypic monitoring of R5-to-X4 transition by motif analysis of human immunodeficiency virus type 1 env V3 loop sequences. J Virol. 2003 Dec;77(24):13376-88. doi: 10.1128/jvi.77.24.13376-13388.2003. PMID: 14645592; PMCID: PMC296044.

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. 2016. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.

J L. Plotrix: a package in the red light district of R. R-News. 2006; 6(4), 8-12.
