#!/bin/bash

#PBS -S /bin/sh
#PBS -N NanOlympicsMod
#PBS -l select=1:ncpus=1:mem=2G
#PBS -M emailAddress
#PBS -m e

source /path/to/activate /path/to/nextflow/environment
cd /path/to/pipeline/folder
nextflow -c pipeline.conf run pipeline.nf -w /path/to/work/folder/
source /path/to/deactivate
