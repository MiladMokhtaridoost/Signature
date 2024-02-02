# Signature Repository

Welcome to the Signature repository! Signature is an innovative machine learning method to comprehensively detect chromosomal interactions, focusing specifically on non-homologous chromosomal contacts (NHCCs). Feel free to explore each folder for a comprehensive understanding of the Signature method. If you have questions, encounter issues, or wish to contribute, please don't hesitate to reach out through the repository's communication channels and Maass lab contact information (https://lab.research.sickkids.ca/maass/). Thank you for your interest in Signature!




The repository is structured as follows:

## 1. Signature Folder
The Signature folder encompasses three key subfolders:

### a. Pre-processing
The Pre-processing folder includes scripts essential for preparing raw Hi-C data for our algorithms. This step is crucial in transforming data into a format suitable for further analysis.

### b. LWLR (Supervised Learning)
The LWLR folder contains the Local Weighted Linear Regression (LWLR) approach, a supervised learning method within Signature. This approach systematically assesses spatial interactions and their extent across diverse cell types. Detailed information and implementation instructions are provided within the LWLR folder.

### c. Community Detection (Unsupervised Learning)
The Community Detection folder houses scripts related to the unsupervised learning approach in Signature, specifically the Community Detection method. This method aims to uncover patterns of chromosomal interactions without prior labeling. Explore this folder for insights into the Community Detection approach and guidance on implementation.

### d. Demo data
This folder contains Signature's input data (cooler interaction frequencies) for two Hi-C datasets. These files can be used as toy data for running a demo of Signature.

## 2. Figures Folder
The Figures folder is dedicated to scripts associated with the figures presented in the "Inter-chromosomal contacts demarcate genome topology along a spatial gradient" paper. These scripts involve a certain level of data pre-processing to generate the figures showcased in our research paper. Refer to this folder for insights into the data processing steps and to reproduce the figures.
