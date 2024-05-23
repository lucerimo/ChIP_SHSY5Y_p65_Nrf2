#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name="QCchip"
#SBATCH --output=test.out.%j


cd /home/lucerimo/chipseq

source activate bioinfo
module load fastqc-0.11.9
module load trimmomatic-0.39

# Crear una variable con nombre filename que lea el archivo que acabamos de crear

filename='samples'

# Crear las variables de entrada para fastqc

DATA=/home/lucerimo/chipseq/RAW
RFQC=/home/lucerimo/chipseq/fastqc_raw
OUTPUT=/home/lucerimo/chipseq/TRIMMED
TFQC=/home/lucerimo/chipseq/fastqc_trimm

# El loop va a leer filename una linea a la vez
while read SAMPLE;
do
    #echo "FastQC ${SAMPLE}"  
    # Crear los nombres de los files
    F=$DATA/${SAMPLE}_1.fq.gz
    R=$DATA/${SAMPLE}_2.fq.gz
    echo "Trimming ${SAMPLE}" 
    FP=$DATA/${SAMPLE}forward_paired.fq.gz
    RP=$DATA/${SAMPLE}reverse_paired.fq.gz
    FU=$DATA/${SAMPLE}forward_unpaired.fq.gz
    RU=$DATA/${SAMPLE}reverse_unpaired.fq.gz
    
# Correr el fastqc para los reads crudos, trimming y QC de los trimmed.
  fastqc -t 5 $F $R -o $RFQC
  trimmomatic PE -threads 20 $F $R $FP $FU $RP $RU ILLUMINACLIP:adapters-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 
  fastqc -t 5 $FP $RP -o $TFQC
  rm $F
  rm $R
  rm $FU
  rm $RU
  mv $FP $OUTPUT
  mv $RP $OUTPUT
done <  $filename 

mail-report luis.riveramontero@ucr.ac.cr 
