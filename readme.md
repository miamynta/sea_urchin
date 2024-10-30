# R scripts and data of sea urchin intestinal bacterial community analysis

This manuscripts contains amplicon raw data and an R script for analyzing data for the following manuscript:

**Sea urchin intestinal bacterial communities depend on seaweed diet and contain nitrogen-fixing symbionts**

by

Mia Bengtsson, Marita Helgesen, Haitao Wang, Stein Fredriksen, Kjell Magnus Norderhaug

Contact for this repository: mia.bengtsson@uni-greifswald.de

## Data files contained in this repository

   * `urchin_map.csv` - contextual data including information on treatments
   * `urchin_seqtab.csv` - ASV table
   * `urchin_tax.csv` - taxonomic classification of ASVs

## Software dependencies

   * R version ...
   * R-packages: vegan, devtools, pairwiseAdonis, phyloseq, DESeq2, ggplot2

## Data analysis documentation

The file urchin_microbiome_analysis.R contains all commands that were used to analyse the sea urchin microbiome data.

## How to obtain a copy of this repository

On your local machine with git installed, execute in a terminal:

```
git clone https://github.com/miamynta/sea_urchin.git
```

This will create a local folder sea_urchin with all contents of this repository.
