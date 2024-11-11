# Signature Repository

Welcome to the Signature repository! Signature is an innovative machine learning method to comprehensively detect chromosomal interactions, focusing specifically on non-homologous chromosomal contacts (NHCCs). Feel free to explore each folder for a comprehensive understanding of the Signature method. If you have questions, encounter issues, or wish to contribute, please don't hesitate to reach out through the repository's communication channels or the [Maass lab website](https://lab.research.sickkids.ca/maass/). Thank you for your interest in Signature!

If you use Signature in your work, please cite:
Mokhtaridoost, M. et al. Signature https://doi.org/10.5281/zenodo.13873973 (2024).

<br/>

The repository is structured as follows:  
<br/>

## 1. Pipeline Folder
The Signature pipeline folder encompasses four key subfolders:

### a. Pre-processing
The Pre-processing folder includes scripts essential for preparing raw Hi-C data for our algorithms. This step is crucial in transforming data into a format suitable for further analysis.

### b. LWPR (Supervised Learning)
The LWPR folder contains the Local Weighted Polynomial Regression (LWPR) approach, a supervised learning method within Signature. This approach systematically assesses spatial interactions and their extent across diverse cell types. Detailed information and implementation instructions are provided within the LWPR folder.

### c. Community Detection (Unsupervised Learning)
The Community Detection folder houses scripts related to the unsupervised learning approach in Signature, specifically the Community Detection method. This method aims to uncover patterns of chromosomal interactions without prior labeling. Explore this folder for insights into the Community Detection approach and guidance on implementation.

### d. Demo
This folder contains Signature's input data (cooler interaction frequencies) for two Hi-C datasets. These files can be used as toy data for running a demo of Signature.  
<br/>

## 2. Figures Folder
The Figures folder is dedicated to scripts associated with the figures presented in the "Inter-chromosomal contacts demarcate genome topology along a spatial gradient" paper. These scripts involve a certain level of data pre-processing to generate the figures showcased in our research paper. Refer to this folder for insights into the data processing steps and to reproduce the figures.
