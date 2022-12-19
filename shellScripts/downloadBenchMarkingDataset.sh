#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH -p cpu
#SBATCH --mem 50G

PATH=$R_HOME/bin:$PATH

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/HG002

wget -co ultra-long-ont.fastq.gz ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V2.3.4_2019-06-26/ultra-long-ont.fastq.gz