# Diploid-GWAS

A workflow for quality control (QC), GWAS, and annotation for diploid species.

Genome-wide association studies require previous quality control and subsequent annotation analyses of genes close to the markers. There are tools to perform these analyses in the R/Bioconductor. However, some of them are no available for non-human data. Additionally, the input of the data may take a long time. This workflow facilitates the input of non-human data with the GeneSeek structure and with a general structure. After this step, there are automatically generates QC, GWAS, and annotation reports. There was used the packages GWASTools, ade4, and SNPRelate.

Code and toy data are in three separate folders: QC and GWAS for data with GeneSeek structure, QC and GWAS for data with general structure, and annotation.
