# PhD_scripts

## Scripts used to generate the work found in my thesis

## Sections
Work is split into two branches, the Drosophila analysis used in chapter 4 of my thesis, and the human analysis used in chapter 5. The sections are all ordered numerically in the order that they must be run (with the exception of functions to source and plotting). These sections are further divided as below:

## Drosophila_Scripts
Unlike in humans (described below) there were only two cell lines used in the Drosophila analysis found in this thesis. As such, these scripts are not split by cell type as in the human directory, but all apply to BG3 and S2.

### 1_Preprocessing
Contains:
1. Scripts used to tile and fill metadata for the BG3 and S2 genomes with metadata read in using the import function from rtracklayer
2. 1 million bin sampling of each genome

### 2_Model_Output_Processing
Contains:
1. Scripts used to read in the .tsv files from the neural netowrk and fuzzy logic model predictions across the entire genome
2. Scripts used to grow the predicted regions

### 3_Enriched_Contacts
Contains:
1. The script used to annotate promoter regions using the .gtf file
2. The scripts used to create GRange objects at each enriched contact that contains a promoter
3. Generation of figure 3A for the paper, and creates a list of enhancers by contact class (Distal Only, Proximal Only, Distal and Proximal, and Neither)
4. Two scripts, one for BG3 and one for S2, these create an object which lists distal and proximal contacts
5. Scripts that create grangelist objects that pair each enhancer with all of its proximal and distal contacts. This is performed in both cell lines for Background, Common, Putative, and STARR-seq only enhancers. A matrix used in future analysis is also created that counts the number of enhancers in each group that do not make contact with any promoters.
6. The expression analysis script that runs on each grangelist created in step 5. This creates the Rda objects that will be used to generate the contact and expression plots 4.10 and 4.11. These objects are also used to generate the plots found in figures S4.1 and S4.2.

### 4_Developmental_vs_housekeeping
Contains:
1. Script used to create the developmental vs housekeeping enhancer dataset
2. Sampling of the developmental vs housekeeping dataset to build models on the Temenos XAI platform

### 5_Arch_Proteins
Contains:
1. The script used to prepare the dataset that Mila created to investigate the additional epigenetic features for uploading to Temenos

### 6_PWM_Enrich
Contains:
- A range of scripts used to query PWMEnrich to create the files used in figure 4.29, as well as additional analysis not included in this thesis

### Functions_to_source
Contains:
- Functions sourced in other scripts

### Plots
Contains:
- Functions to generate each plot found in chapter 4 this thesis

## Human_Scripts
In humans, several of the same analyses were performed in multiple cell types, or with additional unpublished tracks. The following sections detail how these are divided, and what scripts are available for each cell line.

## H9_and_H1
Scripts used in the H9 and H1 analysis. Where different scripts exist for each cell line, the script associated with H9 cells is labelled A and the scripts associated with H1 cells are labelled B.
Contains the following:

### 1_Pre-processing
Contains:
1. Scripts to tile the H9 and H1 genomes into 100bp long GRange objects and populate the metadata columns with metadata read in using the import function from rtracklayer
2. Scripts to sample the H9 and H1 datasets for training on the Temenos XAI platform
3. A script to create a STARR-seq GRanges object for ease of use in future analysis

### 2_Processing_Predicitons
Contains:
1. Scripts to read in the .tsv files containing the XAI predictions for H9 and H1 cells produced by the Temenos platform
2. Scripts to combine enhancers predicted above the threshold of 0.8 within a range of minimum distances for further analysis
3. Scripts to create a full dataset for both cell lines, this includes enhancer annotation (putative, common, STARR-seq only, or neither), annotated region type, sense gene name and antisense gene name. Additional metadata columns not used in the model are removed at this step, and enhancers are reduced with a threshold of 0.8, a minimum.gapwidth of 1000, and a minimum length of 200bp or greater (as described in methods) before enhancer annotation is determined.

### 3_Enriched_Contact_Scripts
Contains:
1. A script to read in Hi-C data processed by Liv and to generate a list of enriched contacts for each enhancer annotation in each cell line
2. Expression processing scripts for each enhancer annotation group in each cell line

### 4_Motif_Enrichment
Contains:
- A script to perform motif enrichment using PWMEnrich in H9 cells

## HeLa-S3
Scripts used in the HeLa-S3. Scripts associated with the strict enhancer definition are labelled A, while those associated withteh lenient enhancer definition are labelled B.
Contains the following:
1. Tiled genome construction of the HeLa-S3 dataset with metadata read in using the import function from rtracklayer
2. Sampling of the two HeLa-S3 datasets to be used in training models on the Temenos XAI platform

## H9_Unpublished_Tracks
Scripts used to add the additional unpublished features in H9 provided by Pradeep.
1. Tiled genome construction of the H9 cell line with metadata read in using the import function from rtracklayer
2. Sampling of the expanded H9 dataset for training models on the Temenos XAI platform

### Functions_to_source
Contains:
- Functions sourced in other scripts

### Plots
Contains:
- Functions to generate each plot found in chapter 5 this thesis
