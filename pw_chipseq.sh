#!/bin/bash

echo ""
echo "[pw_chipseq]"
echo ""


inFileBase=$1 # Basename of input file

# Echo input
INDEX="/proj/dangl_lab/sbiswas/interactome/ref/tair10/indexes/Arabidopsis_thaliana_full_genome_TAIR10_transgene_HopBB1_Thu_Dec_19_18-37-12_EST_2013"
echo "#### INPUT ####"
echo "Input file basename = $1"
echo "Index = $INDEX"
echo ""


infq=$1.fastq
hqfq=$1_HQ.fastq
sam=$1_HQ.sam
ename=$1_HQ

# Echo output
echo "#### OUTPUT ####"
echo "Assumed input file FASTQ = $infq"
echo "Quality filtered FASTQ = $hqfq"
echo "Mapped SAM = $sam"
echo "MACS2 expt. name = $ename"

echo ""


## Quality filter
echo "Running quality filtering ... "
cmd="trimmomatic SE -threads 6 -phred33 $infq $hqfq ILLUMINACLIP:/nas02/apps/trimmomatic-0.32/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:25 TRAILING:20 LEADING:20 MINLEN:36"
echo "EXECUTING: $cmd"
$cmd

## Map
echo "Mapping ..."
cmd="bowtie2 -x /proj/dangl_lab/sbiswas/interactome/ref/tair10/indexes/Arabidopsis_thaliana_full_genome_TAIR10_transgene_HopBB1_Thu_Dec_19_18-37-12_EST_2013 -U $hqfq -S $sam"
echo "EXECUTING: $cmd"
$cmd

## Call peaks
echo "Peak calling ..."
cmd="macs2 callpeak -t $sam -f SAM -g 1.19e8 -n $ename -B"
echo "EXECUTING: $cmd"
$cmd