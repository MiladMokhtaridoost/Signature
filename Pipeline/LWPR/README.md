# Signature: A computational pipeline for determining significant genomic interactions
Signature analyzes Hi-C/Omni-C data to estimate significant interactions for both intra- and inter-chromosomal interactions genome-wide or for a given region of interest. There are three main requirements to run, the pre-processed datasets you want to analyze, the genomic resolution, and the type of analysis. This version of Signature is not prepared to complete trans analyses at resolutions higher than 500KB, or cis analyses at resolutions higher than 50KB. This pipeline is made to run on the Slurm Linux system - the job-submission scripting will need to be modified in all shell scripts to match your High Performance Computing system's guidelines.  
<br/>

## Requirements
### Runs on the Slurm Linux system (modify to your High Performance Computing system guidelines)
Rough guidelines for requirements, based on running 3 datasets, highest resolutions available:
   - Inter-chromosomal analysis (trans, max 500 Kb): ~16Gb memory, ~12hr runtime
   - Intra-chromosomal analysis (cis, max 50 Kb): ~100Gb memory, ~3 days runtime
### Modules
   - bedtools/2.24.0
   - R/4.2.1
   - python/2.7.9
<br/>

## Folder contents
### batch_create_signature_files.sh
   - user-friendly shell script required for running Signature
   - all-in-one script that will generate a configuration file and Slurm schedular to run Signature
   - minimal edits are required (4 lines for _cell info_, 6 lines for _analysis info_)
### batch_create_signature_files_DemoExample.sh
   - filled out example for the toy data from _Signature/Pipeline/Demo_
   - note: this script works for 1-3 cells without modification, for 4+ cells see _Making user-specific modifications_ section below)
### scripts_Signature
   - folder containing all scripts required to run Signature
      1. main shell script with the supervised machine learning (*run_Signature.sh*)
      2. dependency scripts for machine learning (*7 R, 1 Awk, 1 shell*)
<br/>

## Running interaction analysis with Signature
Signature runs with one simple shell script - the *user-friendly scheduler* - which has been generated to submit your job in batches, up to 3 datasets per batch. Larger batches can be submitted depending on your HPC's processing power and memory limits (see *Making user-specific modifications* section). 

### To get started:
1. In https://github.com/MaassLab/Signature, press the ***<>Code*** drop-down button
2. Press ***Download ZIP***
3. Extract ZIP folder (**"Signature-main.zip"**) to your directory of choice 
4. Inside *"Signature-main"*, navigate to the directory *"Pipeline"* then *"LWPR"*
5. Within this directory, you will find everything required to run Signature's LWPR

Note: This pathway is now referred to as "path" in line 19 of **"batch_create_signature_files.sh"**

### Filling in the script:
There are only 2 sections of **"batch_create_signature_files.sh"** that need to be filled in. 

**Cell info section**
- *batch* is any numerical value (making sure that each new batch of cells is assigned a unique batch number)
- Ensure that the dataset names match the folder names in the Cooler output

**Analysis info section**
- *path* is the full pathway to the LWPR folder inside the directory that "Signature-main.zip" was extracted in (*/.../Signature-main/Pipeline/LWPR*)
- *coolerpath* is the full pathway to the folder where the saved cooler data is
- *analysis* can be either **cis** for intra-chromosomal interactions or **trans** for inter-chromosomal interactions
- *trans* can be either **1vsAll** indicating a genome-to-genome comparison (the default choice for trans) or **pairwise** indicating a chromosome-to-chromosome comparison
- *resolution* can be one of the following **1000000**, **500000**, **250000**, **100000**, **50000** (only **1000000** or **500000** for trans)
- *res* is your chosen resolution condensed (i.e. **1MB**, **500KB**, **250KB**, **100KB**, **50KB**)
<br/>
 
## Making user-specific modifications
The *user-friendly scheduler* is made to process signature in batches of 3, but can be edited to process more.

**The following adjustments will need to be made**
1. In the shell script, below line 15, add a new line for each additional cell you want to include (following the same format as line 15) <br/>
   - Example with two additional datasets: <br/>
   **cell4=** <br/>
   **cell5=**
2. In the shell script, below line 42, as you did in step 1, add corresponding lines for each new cell you want to include (following the same format as line 42) <br/>
   - Example with two additional datasets: <br/>
   **echo $cell4 >> $path/batch$batch\_cells.txt** <br/>
   **echo $cell5 >> $path/batch$batch\_cells.txt**
3. String replace "$cell1.$cell2.$cell3" with your new cell-variable string
   - Example using *sed* in command line for a batch of 5 cells: <br/>
   sed -i 's/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.**$cell4.$cell5**/g' batch_create_signature_files.sh
<br/>

## Output
After submitting batch_create_signature_files.sh, a new folder will be generated inside LWPR called **results** (this is where Signature does its processing). All files found in **results** are automatically generated by Signature.

### Results folder content
*/.../Signature-main/Pipeline/LWPR/results* will contain the following folders after run completion:

**schedulers**
   - Directory for Signature's input scripts
   - Folders for each batch containing:
      - Configuration file
      - Slurm schedular
      - Slurm standard output file

**output**
   - Directory for Signature's output files (z-scores, p-values, and q-values)
   - Folders for each batch containing tab-separated text files of interacting regions
      - [cells.analysis].zscoresall.txt $\hspace{12pt}$   # all possible interactions and their z-scores, cell types merged into single file
      - [cells.analysis].zscores.txt    $\hspace{24pt}$   # significant (p < 0.05) interactions and their z-scores, cell types merged into single file
      - [cells.analysis].pvalueAll.txt  $\hspace{12pt}$   # all possible interactions and their p-values, cell types merged into single file
      - [cells.analysis].pvalue.txt     $\hspace{30pt}$   # significant (p < 0.05) interactions and their p-values, cell types merged into single file
      - [cells.analysis].qvalueAll.txt  $\hspace{16pt}$   # all possible interactions and their q-values, cell types merged into single file
      - [cells.analysis].qvalue.txt     $\hspace{30pt}$   # significant (q < 0.05) interactions and their q-values, cell types merged into single file
      - [cell.analysis].pvalue_neg.txt  $\hspace{12pt}$   # significant negative (p < 0.05, z < 0) interactions and their p-values, one per cell type
      - [cell.analysis].pvalue_pos.txt  $\hspace{12pt}$   # significant positive (p < 0.05, z > 0) interactions and their p-values, one per cell type
      - [cell.analysis].qvalue_neg.txt  $\hspace{12pt}$   # significant negative (q < 0.05, z < 0) interactions and their q-values, one per cell type
      - [cell.analysis].qvalue_pos.txt  $\hspace{12pt}$   # significant positive (q < 0.05, z > 0) interactions and their q-values, one per cell type

**Example Astrocyte_Spine.H9hESC_day00_Zhang..trans1vsAll.1MB.zscoresall.txt**

   | ID                                                | Astrocyte_Spine       | H9hESC_day00_Zhang     |
   |:--------------------------------------------------|:----------------------|:-----------------------|
   | Achr1.10000000.11000000Bchr10.0.1000000           | 1.420                 | 1.048                  |
   | Achr1.10000000.11000000Bchr10.100000000.101000000 | 1.752                 | 1.691                  |
   | Achr1.10000000.11000000Bchr10.10000000.11000000   | -1.181                | -0.499                 |

ID format example: Achr1.10000000.11000000Bchr10.0.1000000 

- interaction between chromosome 1 at 10,000,000-11,000,000 bp and chromosmoe 10 at 0-1,000,000 bp
