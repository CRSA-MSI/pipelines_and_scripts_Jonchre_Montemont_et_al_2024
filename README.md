# Pipelines and scripts Jonchere Montemont et al

This repository comprise two folders containig part of the pipelines and scripts used to do the analysis in the article Massive Early Perturbation of Alternative Splicing Partly Due to Microsatellite Instability at U2AF-Binding Polypyrimidic Tract Sites During Initiation of Colorectal Cancer.

The Pipeline folder contains nextflow pipelines to pretreat raw sequencing data :
wes\_map\_sort\_option\_trim.nf : used to do trimming, alignment and sorting of bam files with paired fastq.gz files as input.
wes\_pretreat\_fastq\_UMI\_for\_consensus.nf : used to do the trimming, aligment and deduplication by generating consensus reads for raw paired fastq files with UMI. Compared to the previous pipeline, it needs additional steps to extract UMI, groups reads by aligment position and UMI sequence and then generate consensus reads per group that will be used for the final aligment.

The Scripts folder contains R code :


Firstly, the script "Visualisation.7b" and "0.Plot.t.test.32" provide the list of the significantly deregulated exons (skipping) in MSI cancers compared to normal and MSS, with the statistical analysis. 

Secondly, the script "3.analyse.b.R" provide the list of the candidate microsatellites potentially involved in splicing defects (MS.Candidats.Splicing.txt), i.e. the microsatellites present in the acceptor splicing regulatory sequences, and allows you to obtain the list of the 96 MS which instability is correlated with exon skipping, with the statistical analysis.
This script also generates most of the figures.

Finally, these two scripts "2023.06.14_Splicing_EarlyMSI_TIMSI_splicoter_hg19_v4" & "2024.01.26_Splicing_EarlyMSI_TIMSI_splicoter_hg19_v1_REVIEW" provide the study of the dynamic evolution of microsatellite instability in MMR-deficient crypts, adenoma and adenocarcinoma. It is split in two parts, one for the main analysis and the other for the complementary reviewing work.


