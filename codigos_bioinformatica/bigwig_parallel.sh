#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name="bigwig"
#SBATCH --output=test.out.%j


cd /home/lucerimo/chipseq/ALIGN/BAM

source /opt/BIO/deepTools/bin/activate

# Crear una variable con nombre filename que lea el archivo que acabamos de crear

filename='samples1'

# Crear las variables de entrada para DeepTools

DATA=/home/lucerimo/chipseq/ALIGN/BAM


# El loop va a leer filename una linea a la vez
while read SAMPLE;
do
    echo "Aligning ${SAMPLE}"  
    # Crear los nombres de los files
    BAM=/home/lucerimo/chipseq/ALIGN/BAM/${SAMPLE}.filtered_final.bam
    OUTPUT=/home/lucerimo/chipseq/ALIGN/BAM/${SAMPLE}_RPKM.bw
    OUTPUT2=/home/lucerimo/chipseq/ALIGN/BAM/${SAMPLE}_RPGC.bw

# Correr el el bamCoverage. Normaliza de dos maneras. Va a trabajar en 20 cores al mismo tiempo.
  bamCoverage -bs 20 --normalizeUsing RPKM -b $BAM -o $OUTPUT -p 20
  bamCoverage -bs 20 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -b $BAM -o $OUTPUT2 -p 20

done <  $filename  


mail-report luis.riveramontero@ucr.ac.cr 
