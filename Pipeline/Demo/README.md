# Toy data for running Signature's LWLR and Community Detection (Demo)
This "Demo" folder will act as your output directory from running _cooler_. With _cooler_, all results are generated in a seperate folder for each dataset/cell type analyzed. In this case, 2 datasets would have been processed separetly: H9hESCs and Astrocytes.

## LWLR
Signature's LWLR can process both cis and trans data. The same script can be used for each (with appropriate analysis-specific modifications). In Demo, you will see toy data at the 1Mb resolution for both cis and trans. These files are the input for Signature and are named "[analysis].1000000_iced.sorted.txt"
### Instructions for running demo
1. In https://github.com/MaassLab/Signature, press the _**<>Code**_ drop-down button
2. Press _**Download ZIP**_
3. Extract ZIP folder (_"Signature-main.zip"_) to your directory of choice 
4. Inside _"Signature-main"_, navigate to the directory _"Pipeline"_ then _"LWLR"_
5. Within this directory, you will find everything required to run Signature's LWLR

Notes: 
- This pathway is now referred to as "path" in downstream steps (i.e. line 19 of "batch_create_signature_files.sh")
- 
