#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name="Alignment"
#SBATCH --output=test.out.%j


cd /home/lucerimo/chipseq/

source activate bioinfo
module load bowtie2-2.4.5
module load samtools-1.15
module load bedtools-2.29.1
module load sambamba-1.0.1

# Crear una variable con nombre filename que lea el archivo que acabamos de crear

filename='samples'

# Crear las variables de entrada para bowtiw2

IDX=/home/rebeca.campos/genomes/hg38.fa
DATA=/home/lucerimo/chipseq/TRIMMED


# El loop va a leer filename una linea a la vez
while read SAMPLE;
do
    echo "Aligning ${SAMPLE}"  
    # Crear los nombres de los files
    F=$DATA/${SAMPLE}forward_paired.fq.gz
    R=$DATA/${SAMPLE}reverse_paired.fq.gz
    BAM=/home/lucerimo/chipseq/ALIGN/BAM/${SAMPLE}.bam
    OUTPUT=/home/lucerimo/chipseq/ALIGN/BAM/${SAMPLE}.sam
    LOG=/home/lucerimo/chipseq/ALIGN/LOG/${SAMPLE}.log


# Correr el alineador, convertir el SAM a BAM, hacer el sort y el index del BAM. Va a trabajar en 20 cores al mismo tiempo.
    bowtie2 -p 20 -q -x $IDX -1 $F -2 $R -S $OUTPUT 2> $LOG
    echo "Sorting sam"
    samtools view $OUTPUT -h -S -b -o $BAM 
    rm ./ALIGN/BAM/${SAMPLE}.sam
    samtools sort $BAM -o ./ALIGN/BAM/${SAMPLE}_sorted.bam

done <  $filename  

mail-report luis.riveramontero@ucr.ac.cr 
