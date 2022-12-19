#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem 200G
##SBATCH --nodes 11

PATH=$R_HOME/bin:$PATH

source activate miniStructuralVariantAnalysis 

## Concatenate all fastq files into a single file
cat /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/*.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq

## Run fastQC quality check
fastqc /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/*.fastq --outdir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/QualityControl/FastQC

## Run pycoQC quality check
pycoQC -f /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/sequencing_summary.txt -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/QualityControl/PycoQC/pycoQC.html

## Run MinION_QC quality check
 MinIONQC.R -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/sequencing_summary.txt -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/QualityControl/MinION_QC

## PoreChop Adapter Removal
porechop -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/AdapterRemoval/PoreChop/PoreChop.fastq --discard_middle

## nanofilt Trim Reads and Filtering
NanoFilt -l 500 --headcrop 10 </scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/AdapterRemoval/PoreChop/PoreChop.fastq> /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq

## Minimap & Miniasm Genome Assembly
minimap2 -x ava-ont /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq | gzip -1 > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz

miniasm -f /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa

awk '/^S/{print ">"$2"\n"$3}' /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta

## Flye Genome Assembly
flye --nano-raw /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq --genome-size 1m --out-dir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye

## Shasta Genome Assembly
sed -n '1~4s/^@/>/p;2~4p' /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/all_guppy.fasta

/scratch/prj/sgdp_nanopore/software/shasta-Linux-0.11.1 --input  /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/all_guppy.fasta --config Nanopore-UL-May2022

## Assembly-Stats Check Assembly
assembly-stats /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/AssemblyStats/MinimapAndMiniasm.txt
assembly-stats /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye/assembly.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/AssemblyStats/Flye.txt
assembly-stats /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/all_guppy.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/AssemblyStats/Shasta.txt

## Align to Reference Genome
cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/Dnadiff/MinimapMiniasm
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/Dnadiff/Flye
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye/assembly.fasta

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/Dnadiff/Shasta
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/all_guppy.fasta

## Racon Error Correction
minimap2 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/minimap.racon.paf

racon /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/minimap.racon.paf /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

## Minipolish Error Correction
