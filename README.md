# GPCR_structural

For use with the manuscript:

Title:
Oncogenic rewiring of Class A GPCR signaling pathways within human malignancy using recurrent mutations to cognate motif residues
Authors:
Jonathan Gallion*, Sara Wright&, Theodore Wensel$,&, Olivier Lichtarge$
 
Author affiliation:
$Program in Quantitative and Computational Biosciences
&Department of Biochemistry, Baylor College of Medicine, Houston, TX 77030. 

Corresponding authors:
Olivier Lichtarge
One Baylor Plaza, Room T921, MS BCM225
Houston, TX 77030
Phone: (713) 798-5646




# Purpose
The intent of this repository is to make available code to recreate our analysis of modeling mutations from all Class A GPCRS and identify structural hotspots across the entire class.

# Use
The code pipeline is as follows (note: many of the file locations and output files are hard coded to use a specific file structure, you will need to mimic this structure or alter locations) I have provided the smaller input files required to run (e.g. GPCR protein alignment used, naming conversion file, fasta sequences of GPCRS, etc).:

1)Identifying_Frequently_Mutated_Positions.ipynb     This jupyter notebook parses over the raw input data from the TCGA and maps all mutations within the Class A GPCRS specified using the protein alignment file (.fa). Output is a txt file listing all the mutations and their aligned position.
2)Make_pymolfile_for_sphere_sizes.ipynb    This jupyter notebook takes the file from above and maps all mutations onto a a template structure (In this case Im using the inactive structure of B2AR (2RH1)). The output is a .pml file. This can be dragged into pymol in order to format the structure.
3)RandomExpectations_startingfromcodon.py This python script performs a random simulation in to mimic random mutation rates with the assumption that mutations follow the transition tranversion frequencie for that cancer type and will therefore impact different codons differently.
4)ETEAcompare.py   This python script compares ET and EA distributions of the hyper, hypomutated positions compared to random expectation.
5)AAtransfrequencies_simpleplot.py   This python script plots the amino acid substitutions at the hyper mutated positions to identify the primary changes.
6) Gprot_enrichment.py   Tests the hypothesis that mutations may be enriched within specified G-protein pathways.



