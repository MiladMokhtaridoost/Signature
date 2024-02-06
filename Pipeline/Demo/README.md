# Toy data for running Signature's LWLR and Community Detection (Demo)
This "Demo" folder will act as your output directory from running _Cooler_. With _Cooler_, all results are generated in a separate folder for each dataset/cell type analyzed. In this case, 2 cell types would have been processed separately: H9hESCs and Astrocytes.  
<br/>

## LWLR
Signature's LWLR can process both cis and trans data. The same script can be used for each (with appropriate analysis-specific modifications). The Demo directory contains toy data at the 1MB resolution for both cis and trans (input for Signature). These files are named "[analysis].1000000_iced.sorted.txt".

### Instructions for running demo
1. In https://github.com/MaassLab/Signature, press the <>Code drop-down button
2. Press Download ZIP
3. Extract ZIP folder (_"Signature-main.zip"_) to your directory of choice 
4. In the command line, navigate to *"Signature-main"*, then the directory *"Pipeline"*, then *"LWLR"*
5. Fill in the script *"batch_create_signature_files.sh"* using the following information
   -  this is the **first** batch of cells, which is processing 2 cell types with the dataset names **"Astrocyte_Spine"** and **"H9hESC_day00_Zhang"**
   -  the pathway is **"*/your/extracted/location*/Signature-main/Pipeline/LWLR"**
   -  the cooler output data can be found in **"*/your/extracted/location*/Signature-main/Pipeline/Demo"**
   -  the data will be processed in **"trans"** for a **"1vsAll"** approach at a resolution of **"1000000"** (i.e. **"1MB"**)

Note: the above information from step 5 can be seen filled-in in an example script (/.../Signature-main/Pipeline/LWLR/batch_create_signature_files_DemoExample.sh)  
<br/>

After Signature is done running the demo data, the results can be found in "*/your/extracted/location*/Signature-main/Pipeline/LWLR/results/output" and will look something like this:

**Example Astrocyte_Spine.H9hESC_day00_Zhang..trans1vsAll.1MB.zscoresall.txt**

   | ID                                                | Astrocyte_Spine       | H9hESC_day00_Zhang     |
   |:--------------------------------------------------|:----------------------|:-----------------------|
   | Achr1.10000000.11000000Bchr10.0.1000000           | 1.420                 | 1.048                  |
   | Achr1.10000000.11000000Bchr10.100000000.101000000 | 1.752                 | 1.691                  |
   | Achr1.10000000.11000000Bchr10.10000000.11000000   | -1.181                | -0.499                 |

Signature ran this toy data in about 10 mins utilizing 0.52 GB of memory (note: toy data contains ~4-5% of the original data)

<br/>
   
## Community Detection
We established Signature's Community Detection in a way that process cis and trans data alltogether. However, users can use it to analyze only cis or only trans. The Demo directory contains toy data at the 1MB resolution (input for Community Detection). These files are named "[cell].pairs.res1000000.cool.txt" and each are located in their corresponding cell type folder.

### Instructions for running demo
1. In https://github.com/MaassLab/Signature, press the <>Code drop-down button
2. Press Download ZIP
3. Extract ZIP folder (_"Signature-main.zip"_) to your directory of choice 
4. Navigate to *"Signature-main"*, then the directory *"Pipeline"*, then *"Community_Detection"*
5. First, run data_processing.R script without any changes.
   -  This can be run on Demo datasets very fast using any computer. We suggest to use HPC for full datasets running on different datasets in parallel.
   -  The output will be a re-formatted version of cooler data for each dataset, including three columns. First two columns represent the ID of binA and binB in [chr]_[start_of_bin(MB)] format and the third column shows the Hi_C weight. For instance, ID 6_41 refers to chr6:41000000.42000000.  
   -  The outputs will be saved in their corresponding folder (e.g.       **Demo/Astrocyte_Spine/Astrocyte_Spine_network.txt**)
6. Secondly, run merge_networks.R to merge the Hi-C interaction weights of every single bin (cis and trans) of all datasets by taking the average of them.
   -  This can be run on Demo datasets very fast using any computer. Normal computers can handle this code for full collection of datasets (62 datasets), however, for faster performance using HPC is suggested. 
   -  The output will be in a same format as the previous step, but merged all datasets in one file.
   -  The output will be saved in the Demo folder (**Demo/average_1MB_network.txt**)
7. The final step is applying Community Detectipn on the merged Hi-C data. For more information, including requierements please read README file in Community_Detection folder.
   -  This can be run on Demo datasets and full datasets very fast using any computer and HPC. 
   -  The output will be a two-column file, first column shows the ID of a bin and the second column shows the assigned community number. The parameters set in a way that result with full datasets will be 46 communities. However, the Demo file will be resulted in less than 46 cpmmunities.
   -  The final output will be saved in the Demo folder (**Demo/average_netwrok_comms.txt**)
