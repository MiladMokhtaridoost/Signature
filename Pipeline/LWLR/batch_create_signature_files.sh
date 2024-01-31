#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=5G
#SBATCH -t 00:15:00
#SBATCH --output=%x.e%j



#----------------EDIT CELL INFO------------------#
batch=		# number  ex.  [4]

cell1=		# cell_type  ex. [Astrocyte_spine]
cell2=		# cell_type  ex. [Astrocyte_cerebellum]
cell3=		# cell_type  ex. [Aorta_Leung]
#------------------------------------------------#

#---------------EDIT ANALYSIS INFO---------------#
path=			# /pathway/to/folder/this/script/is/in (/.../Signature-main/Pipeline/LWLR)
coolerpath=		# /pathway/to/cooler/output

analysis=      	# either [cis] or [trans]
trans=	      	# if trans selected above: either [1vsAll] or [pairwise] (leave blank if cis was selected above)
resolution=	    # resolution: [1000000], [50000], etc.
res=        	# resolution in human readable format; [1MB], [50KB], etc.
#------------------------------------------------#






#-----------------DONT EDIT----------------------#
scr=$path/scripts_Signature.v29
path=$path/results
mkdir $path

# LIST OF CELLS

echo $cell1 > $path/batch$batch\_cells.txt
echo $cell2 >> $path/batch$batch\_cells.txt
echo $cell3 >> $path/batch$batch\_cells.txt

mkdir -p $path/schedulers/$cell1.$cell2.$cell3
mkdir -p $path/output/$cell1.$cell2.$cell3
echo "Pathway to configuration and scheduler = $path/schedulers/$cell1.$cell2.$cell3"



# CONFIGURATION FILE

true > $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt

echo type.class=$analysis$trans >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo specie=human >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo cell=$path/batch$batch\_cells.txt >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo res=1000000 >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo name=$cell1.$cell2.$cell3.$analysis$trans.$res >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo out=$path/output/$cell1.$cell2.$cell3 >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo pval=0.05 >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo dir=$coolerpath >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo in=$scr >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt

echo "CONFIGURATION FILE DONE"


# SCHEDULER FILE

true > $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo '#!/bin/bash' >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH -N 1" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH -c 20" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH --mem=100G" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH --tmp=300G" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH -p general" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH -t 72:00:00" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH --output=%x.e%j" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo -e "\n" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo "module load R/4.2.1" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "toSign=$scr" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "configuration=$path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo '$toSign/Signature.v29.sh $(cat $configuration)' >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo "SCHEDULER FILE DONE"

echo "SUBMITTING SCHEDULER ..."
cd $path/schedulers/$cell1.$cell2.$cell3
sbatch $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

#------------------------------------------------#
