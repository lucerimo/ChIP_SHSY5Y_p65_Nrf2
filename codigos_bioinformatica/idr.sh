#!/bin/bash
#SBATCH --partition=serial
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="IDR"
#SBATCH --output=test.out.%j


cd /home/lucerimo/chipseq/ALIGN/BAM/IDR

source activate bioinfo
module load idr-2.0.3

#Calculo del irreproducibility discovery rate para cada combinacion pareada.

echo "CP"
idr --samples CP_p65_1_peaks_sorted.narrowPeak CP_p65_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_p65_12-idr.narrowPeak --plot --log-output-file CP_p65_12.idr.log --verbose --idr-threshold 0.05
idr --samples CP_p65_1_peaks_sorted.narrowPeak CP_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_p65_13-idr.narrowPeak --plot --log-output-file CP_p65_13.idr.log --verbose --idr-threshold 0.05
idr --samples CP_p65_1_peaks_sorted.narrowPeak CP_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_p65_14-idr.narrowPeak --plot --log-output-file CP_p65_14.idr.log --verbose --idr-threshold 0.05
idr --samples CP_p65_2_peaks_sorted.narrowPeak CP_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_p65_23-idr.narrowPeak --plot --log-output-file CP_p65_23.idr.log --verbose --idr-threshold 0.05
idr --samples CP_p65_2_peaks_sorted.narrowPeak CP_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_p65_24-idr.narrowPeak --plot --log-output-file CP_p65_24.idr.log --verbose --idr-threshold 0.05
idr --samples CP_p65_3_peaks_sorted.narrowPeak CP_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_p65_34-idr.narrowPeak --plot --log-output-file CP_p65_34.idr.log --verbose --idr-threshold 0.05

idr --samples CP_RNA_1_peaks_sorted.narrowPeak CP_RNA_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_RNA-idr_12.narrowPeak --plot --log-output-file CP_RNA_12.idr.log --verbose --idr-threshold 0.05
idr --samples CP_RNA_1_peaks_sorted.narrowPeak CP_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_RNA-idr_13.narrowPeak --plot --log-output-file CP_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples CP_RNA_2_peaks_sorted.narrowPeak CP_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_RNA-idr_23.narrowPeak --plot --log-output-file CP_RNA_23.idr.log --verbose --idr-threshold 0.05

idr --samples CP_Nrf2_1_peaks_sorted.narrowPeak CP_Nrf2_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_Nrf2-idr_12.narrowPeak --plot --log-output-file CP_Nrf2_12.idr.log --verbose --idr-threshold 0.05
idr --samples CP_Nrf2_1_peaks_sorted.narrowPeak CP_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_Nrf2-idr_13.narrowPeak --plot --log-output-file CP_Nrf2_13.idr.log --verbose --idr-threshold 0.05
idr --samples CP_Nrf2_1_peaks_sorted.narrowPeak CP_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_Nrf2-idr_14.narrowPeak --plot --log-output-file CP_Nrf2_14.idr.log --verbose --idr-threshold 0.05
idr --samples CP_Nrf2_2_peaks_sorted.narrowPeak CP_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_Nrf2-idr_23.narrowPeak --plot --log-output-file CP_Nrf2_23.idr.log --verbose --idr-threshold 0.05
idr --samples CP_Nrf2_2_peaks_sorted.narrowPeak CP_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_Nrf2-idr_24.narrowPeak --plot --log-output-file CP_Nrf2_24.idr.log --verbose --idr-threshold 0.05
idr --samples CP_Nrf2_3_peaks_sorted.narrowPeak CP_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CP_Nrf2-idr_34.narrowPeak --plot --log-output-file CP_Nrf2_34.idr.log --verbose --idr-threshold 0.05


echo "1422"
idr --samples 1422_p65_1_peaks_sorted.narrowPeak 1422_p65_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_p65_12-idr.narrowPeak --plot --log-output-file 1422_p65_12.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_p65_1_peaks_sorted.narrowPeak 1422_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_p65_13-idr.narrowPeak --plot --log-output-file 1422_p65_13.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_p65_1_peaks_sorted.narrowPeak 1422_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_p65_14-idr.narrowPeak --plot --log-output-file 1422_p65_14.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_p65_2_peaks_sorted.narrowPeak 1422_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_p65_23-idr.narrowPeak --plot --log-output-file 1422_p65_23.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_p65_2_peaks_sorted.narrowPeak 1422_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_p65_24-idr.narrowPeak --plot --log-output-file 1422_p65_24.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_p65_3_peaks_sorted.narrowPeak 1422_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_p65_34-idr.narrowPeak --plot --log-output-file 1422_p65_34.idr.log --verbose --idr-threshold 0.05

idr --samples 1422_RNA_1_peaks_sorted.narrowPeak 1422_RNA_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_RNA-idr_12.narrowPeak --plot --log-output-file 1422_RNA_12.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_RNA_1_peaks_sorted.narrowPeak 1422_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_RNA-idr_13.narrowPeak --plot --log-output-file 1422_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_RNA_1_peaks_sorted.narrowPeak 1422_RNA_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_RNA-idr_13.narrowPeak --plot --log-output-file 1422_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_RNA_2_peaks_sorted.narrowPeak 1422_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_RNA-idr_23.narrowPeak --plot --log-output-file 1422_RNA_23.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_RNA_2_peaks_sorted.narrowPeak 1422_RNA_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_RNA-idr_13.narrowPeak --plot --log-output-file 1422_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_RNA_3_peaks_sorted.narrowPeak 1422_RNA_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_RNA-idr_13.narrowPeak --plot --log-output-file 1422_RNA_13.idr.log --verbose --idr-threshold 0.05

idr --samples 1422_Nrf2_1_peaks_sorted.narrowPeak 1422_Nrf2_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_Nrf2-idr_12.narrowPeak --plot --log-output-file 1422_Nrf2_12.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_Nrf2_1_peaks_sorted.narrowPeak 1422_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_Nrf2-idr_13.narrowPeak --plot --log-output-file 1422_Nrf2_13.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_Nrf2_1_peaks_sorted.narrowPeak 1422_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_Nrf2-idr_14.narrowPeak --plot --log-output-file 1422_Nrf2_14.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_Nrf2_2_peaks_sorted.narrowPeak 1422_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_Nrf2-idr_23.narrowPeak --plot --log-output-file 1422_Nrf2_23.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_Nrf2_2_peaks_sorted.narrowPeak 1422_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_Nrf2-idr_24.narrowPeak --plot --log-output-file 1422_Nrf2_24.idr.log --verbose --idr-threshold 0.05
idr --samples 1422_Nrf2_3_peaks_sorted.narrowPeak 1422_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file 1422_Nrf2-idr_34.narrowPeak --plot --log-output-file 1422_Nrf2_34.idr.log --verbose --idr-threshold 0.05


echo "Resv"
idr --samples Resv_p65_1_peaks_sorted.narrowPeak Resv_p65_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_p65_12-idr.narrowPeak --plot --log-output-file Resv_p65_12.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_p65_1_peaks_sorted.narrowPeak Resv_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_p65_13-idr.narrowPeak --plot --log-output-file Resv_p65_13.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_p65_1_peaks_sorted.narrowPeak Resv_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_p65_14-idr.narrowPeak --plot --log-output-file Resv_p65_14.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_p65_2_peaks_sorted.narrowPeak Resv_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_p65_23-idr.narrowPeak --plot --log-output-file Resv_p65_23.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_p65_2_peaks_sorted.narrowPeak Resv_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_p65_24-idr.narrowPeak --plot --log-output-file Resv_p65_24.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_p65_3_peaks_sorted.narrowPeak Resv_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_p65_34-idr.narrowPeak --plot --log-output-file Resv_p65_34.idr.log --verbose --idr-threshold 0.05

idr --samples Resv_RNA_1_peaks_sorted.narrowPeak Resv_RNA_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_RNA-idr_12.narrowPeak --plot --log-output-file Resv_RNA_12.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_RNA_1_peaks_sorted.narrowPeak Resv_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_RNA-idr_13.narrowPeak --plot --log-output-file Resv_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_RNA_1_peaks_sorted.narrowPeak Resv_RNA_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_RNA-idr_13.narrowPeak --plot --log-output-file Resv_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_RNA_2_peaks_sorted.narrowPeak Resv_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_RNA-idr_13.narrowPeak --plot --log-output-file Resv_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_RNA_2_peaks_sorted.narrowPeak Resv_RNA_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_RNA-idr_13.narrowPeak --plot --log-output-file Resv_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_RNA_3_peaks_sorted.narrowPeak Resv_RNA_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_RNA-idr_23.narrowPeak --plot --log-output-file Resv_RNA_23.idr.log --verbose --idr-threshold 0.05

idr --samples Resv_Nrf2_1_peaks_sorted.narrowPeak Resv_Nrf2_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_Nrf2-idr_12.narrowPeak --plot --log-output-file Resv_Nrf2_12.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_Nrf2_1_peaks_sorted.narrowPeak Resv_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_Nrf2-idr_13.narrowPeak --plot --log-output-file Resv_Nrf2_13.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_Nrf2_1_peaks_sorted.narrowPeak Resv_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_Nrf2-idr_14.narrowPeak --plot --log-output-file Resv_Nrf2_14.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_Nrf2_2_peaks_sorted.narrowPeak Resv_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_Nrf2-idr_23.narrowPeak --plot --log-output-file Resv_Nrf2_23.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_Nrf2_2_peaks_sorted.narrowPeak Resv_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_Nrf2-idr_24.narrowPeak --plot --log-output-file Resv_Nrf2_24.idr.log --verbose --idr-threshold 0.05
idr --samples Resv_Nrf2_3_peaks_sorted.narrowPeak Resv_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file Resv_Nrf2-idr_34.narrowPeak --plot --log-output-file Resv_Nrf2_34.idr.log --verbose --idr-threshold 0.05


echo "CN"
idr --samples CN_p65_1_peaks_sorted.narrowPeak CN_p65_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_p65_12-idr.narrowPeak --plot --log-output-file CN_p65_12.idr.log --verbose --idr-threshold 0.05
idr --samples CN_p65_1_peaks_sorted.narrowPeak CN_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_p65_13-idr.narrowPeak --plot --log-output-file CN_p65_13.idr.log --verbose --idr-threshold 0.05
idr --samples CN_p65_1_peaks_sorted.narrowPeak CN_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_p65_14-idr.narrowPeak --plot --log-output-file CN_p65_14.idr.log --verbose --idr-threshold 0.05
idr --samples CN_p65_2_peaks_sorted.narrowPeak CN_p65_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_p65_23-idr.narrowPeak --plot --log-output-file CN_p65_23.idr.log --verbose --idr-threshold 0.05
idr --samples CN_p65_2_peaks_sorted.narrowPeak CN_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_p65_24-idr.narrowPeak --plot --log-output-file CN_p65_24.idr.log --verbose --idr-threshold 0.05
idr --samples CN_p65_3_peaks_sorted.narrowPeak CN_p65_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_p65_34-idr.narrowPeak --plot --log-output-file CN_p65_34.idr.log --verbose --idr-threshold 0.05

idr --samples CN_RNA_1_peaks_sorted.narrowPeak CN_RNA_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_RNA-idr_12.narrowPeak --plot --log-output-file CN_RNA_12.idr.log --verbose --idr-threshold 0.05
idr --samples CN_RNA_1_peaks_sorted.narrowPeak CN_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_RNA-idr_13.narrowPeak --plot --log-output-file CN_RNA_13.idr.log --verbose --idr-threshold 0.05
idr --samples CN_RNA_2_peaks_sorted.narrowPeak CN_RNA_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_RNA-idr_23.narrowPeak --plot --log-output-file CN_RNA_23.idr.log --verbose --idr-threshold 0.05

idr --samples CN_Nrf2_1_peaks_sorted.narrowPeak CN_Nrf2_2_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_Nrf2-idr_12.narrowPeak --plot --log-output-file CN_Nrf2_12.idr.log --verbose --idr-threshold 0.05
idr --samples CN_Nrf2_1_peaks_sorted.narrowPeak CN_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_Nrf2-idr_13.narrowPeak --plot --log-output-file CN_Nrf2_13.idr.log --verbose --idr-threshold 0.05
idr --samples CN_Nrf2_1_peaks_sorted.narrowPeak CN_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_Nrf2-idr_14.narrowPeak --plot --log-output-file CN_Nrf2_14.idr.log --verbose --idr-threshold 0.05
idr --samples CN_Nrf2_2_peaks_sorted.narrowPeak CN_Nrf2_3_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_Nrf2-idr_23.narrowPeak --plot --log-output-file CN_Nrf2_23.idr.log --verbose --idr-threshold 0.05
idr --samples CN_Nrf2_2_peaks_sorted.narrowPeak CN_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_Nrf2-idr_24.narrowPeak --plot --log-output-file CN_Nrf2_24.idr.log --verbose --idr-threshold 0.05
idr --samples CN_Nrf2_3_peaks_sorted.narrowPeak CN_Nrf2_4_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file CN_Nrf2-idr_34.narrowPeak --plot --log-output-file CN_Nrf2_34.idr.log --verbose --idr-threshold 0.05

#reacomodo de los outputs en carpeta respectiva
mv *.png ./png
mv *.log ./log
wc -l *-idr_**.narrowPeak

done <  $filename  

mail-report luis.riveramontero@ucr.ac.cr 
