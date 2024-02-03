# Toy data for running Signature's LWLR and Community Detection (Demo)
This "Demo" folder will act as your output directory from running _Cooler_. With _Cooler_, all results are generated in a seperate folder for each dataset/cell type analyzed. In this case, 2 cell types would have been processed separetly: H9hESCs and Astrocytes.  
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

Note: the above information from step 5 can be seen filled-in in an exmaple script (/.../Signature-main/Pipeline/LWLR/batch_create_signature_files_DemoExample.sh)  
<br/>
   
## Community Detection
MILAD FILL SOMETHING IN HERE. The Demo directory contains toy data at the 1MB resolution (input for Community Detection). These files are named "[cell].pairs.res1000000.cool.txt".

### Instructions for running demo
1. In https://github.com/MaassLab/Signature, press the <>Code drop-down button
2. Press Download ZIP
3. Extract ZIP folder (_"Signature-main.zip"_) to your directory of choice 
4. In the command line, navigate to *"Signature-main"*, then the directory *"Pipeline"*, then *"Community_Detection"*
5. Fill in the script ....?
6. continue here Milad
