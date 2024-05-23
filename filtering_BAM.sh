#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name="Filtering_BAM"
#SBATCH --output=test.out.%j


cd /home/lucerimo/chipseq/ALIGN/BAM

source activate bioinfo
module load bedtools-2.29.1
module load sambamba-1.0.1

# Crear una variable con nombre filename que lea el archivo que acabamos de crear

filename='samples'

# Crear las variables de entrada para bedtools

IDX=/home/lucerimo/chipseq/ALIGN/hg38-blacklist.v2.bed


# El loop va a leer filename una linea a la vez
while read SAMPLE;
do
    echo "Filtering ${SAMPLE}"  
    # Crear los nombres de los files
    BAM=${SAMPLE}_sorted.bam
    FILTERED=${SAMPLE}_filtered.bam
    FINAL=${SAMPLE}.filtered_final.bam


# Filtrar secuencias no mapeadas y en las regiones pertenecientes al blacklist.
    sambamba view -h -t 20 -f bam -F "[XS] == null and not unmapped and not duplicate" $BAM > $FILTERED
    bedtools intersect -v -abam $FILTERED -b $IDX > $FINAL

done <  $filename  



mail-report luis.riveramontero@ucr.ac.cr 
