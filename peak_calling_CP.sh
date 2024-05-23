#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name="Peak_call"
#SBATCH --output=test.out.%j


cd /home/lucerimo/chipseq/ALIGN/BAM

source /opt/BIO/Macs/bin/activate


# Corre MACS para identificar los picos usando las replicas de IgG por condicion biologica.
  echo "NRF2"
    macs3 callpeak -t C58.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_Nrf2_1 --outdir macs3 2> macs3/CP_Nrf2_1_macs3.log
    macs3 callpeak -t C62.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_Nrf2_2 --outdir macs3 2> macs3/CP_Nrf2_2_macs3.log
    macs3 callpeak -t C90.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_Nrf2_3 --outdir macs3 2> macs3/CP_Nrf2_3_macs3.log
    macs3 callpeak -t C94.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_Nrf2_4 --outdir macs3 2> macs3/CP_Nrf2_4_macs3.log
    
    # Sorting de los picos por p(p_value)
    cat CP_Nrf2_1_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_Nrf2_1_sorted.narrowPeak
    cat CP_Nrf2_2_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_Nrf2_2_sorted.narrowPeak
    cat CP_Nrf2_3_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_Nrf2_3_sorted.narrowPeak
    cat CP_Nrf2_4_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_Nrf2_4_sorted.narrowPeak        
    
    
  echo "p65"    

    macs3 callpeak -t C50.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_p65_1 --outdir macs3 2> macs3/CP_p65_1_macs3.log
    macs3 callpeak -t C54.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_p65_2 --outdir macs3 2> macs3/CP_p65_2_macs3.log
    macs3 callpeak -t C82.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_p65_3 --outdir macs3 2> macs3/CP_p65_3_macs3.log
    macs3 callpeak -t C86.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_p65_4 --outdir macs3 2> macs3/CP_p65_4_macs3.log

    # Sorting de los picos por p(p_value) 
    cat CP_p65_1_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_p65_1_sorted.narrowPeak
    cat CP_p65_2_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_p65_2_sorted.narrowPeak
    cat CP_p65_3_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_p65_3_sorted.narrowPeak
    cat CP_p65_4_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_p65_4_sorted.narrowPeak
    
  echo "RNApolII" 
    
    macs3 callpeak -t C38.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_RNA_1 --outdir macs3 2> macs3/CP_RNA_1_macs3.log
    macs3 callpeak -t C66.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_RNA_2 --outdir macs3 2> macs3/CP_RNA_2_macs3.log
    macs3 callpeak -t C70.filtered_final.bam -c C46.filtered_final.bam C74.filtered_final.bam C78.filtered_final.bam  -f BAM -g hs -n CP_RNA_3 --outdir macs3 2> macs3/CP_RNA_3_macs3.log

    # Sorting de los picos por p(p_value)
    cat CP_RNA_1_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_RNA_1_sorted.narrowPeak
    cat CP_RNA_2_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_RNA_2_sorted.narrowPeak
    cat CP_RNA_3_peaks.narrowPeak|sort -k8,8nr|head -n 100000 >CP_RNA_3_sorted.narrowPeak
    
done <  $filename  

mail-report luis.riveramontero@ucr.ac.cr 