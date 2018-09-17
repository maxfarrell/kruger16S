# Analyses for "Bacterial diversity in the waterholes of the Kruger National Park: an eDNA metabarcoding approach"

[10.1139/gen-2018-0064](10.1139/gen-2018-0064)

## Scripts
1. DADA2 pipeline for quality filtering, denoising, chimera removal, ASV calling, and taxonomy assignment (dada2_pipeline.R)
2. Build phylogeny (build_phylo.R) - R script with system calls to QIIME and FastTree
3. Community and biodiversity analyses (community_analysis.R)

### Required data

The community analyses can be reproduced with the data files included here. To reproduce the DADA2 pipeline and phylogeny building the raw sequence reads are archived in the NCBI Sequence Read Archive: 

* BioProject PRJNA490450 
* Accession numbers SRR7822814 to SRR7822901

In addition you will need the [DADA2 formatted training fasta files for SILVA version 128](https://zenodo.org/record/801832#.W5_UPPMbCEI) and the [SILVA 128 release for QIIME](https://www.arb-silva.de/download/archive/qiime/).