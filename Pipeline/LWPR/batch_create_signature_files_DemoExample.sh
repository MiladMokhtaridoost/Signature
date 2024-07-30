#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=5G
#SBATCH -t 00:15:00
#SBATCH --output=%x.e%j



#----------------EDIT CELL INFO------------------#
batch=1

cell1=H9hESC_day00_Zhang
cell2=Astrocyte_Spine
cell3=
#------------------------------------------------#

#---------------EDIT ANALYSIS INFO---------------#
path=/your/hpc/directory/Signature-main/Pipeline/LWLR
coolerpath=/your/hpc/directory/Signature-main/Pipeline/Demo

analysis=trans
trans=1vsAll
resolution=1000000
res=1MB
#------------------------------------------------#

