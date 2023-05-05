# NanOlympicsMod

This repository contains the code developed for benchmarking m6A-detection tools based on Nanopore direct RNA sequencing.

## Getting started

**Prerequisites**

* [Nextflow](https://nf-co.re/usage/installation)
* [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)                                                                                  
                                                                                   
**Installation**

```
git clone https://github.com/mfurla/NanOlympicsMod.git
cd NanOlympicsMod
chmod 755 *
```

## Usage

The NanOlympicsMod pipeline requires you to open pipeline.conf configuration file and set the desired options.
```
   Usage:
   nextflow -c pipeline.conf run pipeline.nf --samples="/path/to/samples.txt" --resultsDir="/path/to/resultsDir" 
   Mandatory arguments which may be specified in the pipeline.conf file

--samples                                                Path to the tab-separated sample file including sample name, condition and path to base-called fast5 folder
--test_condition                                         Condition that we are interested to profile (e. g. 'WT')
--resultsDir                                             Path to a folder where to store results
--fast5_slot                                             FAST5 slot containing the basecalled bases
--fast5_slot_id                                          FAST5 slot containing the basecalled bases (redundant)
--tombo_slot                                             FAST5 slot containing the resquiggled data
--tombo_subslot                                          FAST5 slot containing the resquiggled data
--transcriptome_fasta                                    Path to the transcriptome fasta file
--transcriptome_fai                                      Path to the transcriptome fasta index file
--genome_fasta                                           Path to the genome fasta file
--genome_fai                                             Path to the genome fasta index file
--genes2transcripts                                      Path to gene-to-transcripts file for Nanom6A
--transcriptomebed                                       Path to transcripts bed12 file
--genesbed                                               Path to genes bed file
--gtf                                                    Path to genome annotation gtf file
--nanom6AP                                               nanom6A probability thresholds for PR curve plotting
--yanocompFDR                                            yanocomp FDR threshold
--differrFDR                                             differr FDR threshold
--drummerPval                                            drummer Pvalue threshold
--epinanoErrorSumErr                                     epinanoError threshold sum of errors
--epinanoErrorResiduals                                  epinanoError threshold residuals
--postprocessingScript                                   Path to postprocessing R script
--statisticalAnalysis                                    Path to statistical_analysis R script
--binLength                                              Size of windows for genome binning
--threshold                                              Set of thresholds to use for the filtering of m6A sites (choose between 'default' and 'relaxed') 
--peaksfile                                              Path to bed file with set of m6A gold-standard peaks
```
