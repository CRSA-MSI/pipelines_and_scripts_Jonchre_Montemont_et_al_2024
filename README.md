# Pipelines and scripts Jonchere Montemont et al

This repository comprise two folders containig part of the pipelines and scripts used to do the analysis in the article Massive Early Perturbation of Alternative Splicing Partly Due to Microsatellite Instability at U2AF-Binding Polypyrimidic Tract Sites During Initiation of Colorectal Cancer.

The Pipeline folder contains nextflow pipelines to pretreat raw sequencing data :
wes\_map\_sort\_option\_trim.nf : used to do trimming, alignment and sorting of bam files with paired fastq.gz files as input.
wes\_pretreat\_fastq\_UMI\_for\_consensus.nf : used to do the trimming, aligment and deduplication by generating consensus reads for raw paired fastq files with UMI. Compared to the previous pipeline, it needs additional steps to extract UMI, groups reads by aligment position and UMI sequence and then generate consensus reads per group that will be used for the final aligment.

The Scripts folder contains R code :
