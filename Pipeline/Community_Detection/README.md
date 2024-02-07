# Estimating spatial genome topology with Community Detection

We used Community Detection (CD) to explore spatial genome topology. CD generates clusters of cis and trans interaction frequencies, determines multi-way interactions, and visualizes the CD-derived results in genome topology maps.

#### Requirements
- python/3.8.0
  
  - Pycombo (https://pypi.org/project/pycombo/) 

- R/4.2.1


#### Input Data

To perform CD, we utilized the interaction frequency (including cis and trans) from cooler of our compendium of 62 individual Hi-C datasets (as detailed in the paper) for each interaction.

#### Utilizing Combo Algorithm

We implemented the Combo algorithm, which is available through the pycombo Python package (https://pycombo.readthedocs.io/en/latest/), as a non-overlapping CD algorithm.


#### Output
After the completion of the Community Detection step, each bin will be assigned a community group number. The parameters are tuned to generate 46 distinct communities. Here is an example of output format:

   | ID name    | Community ID       | 
   |:--------------------------------------------------|:----------------------|
   | 1_3  | 1                 | 
   | 1_176  | 23                 |
   | 3_40  | 1                | 
   | 21_32  | 14                | 

For example, ID name 1_176 refers to the genomic region chr1:176000000-177000000. According to this example output 1_3 and 3_40 are in the same community.

## Scripts

### data_processing.R

Data Preparation: In this phase, we prepared the data by converting the interaction frequencies (from cooler) for each interaction (cis and trans) across all datasets into a network format suitable for subsequent CD analysis. Please refer to the Demo folder for more information (including available input toy data). Here is an example of output format:




### merge_netwroks.R
In this step, we are computing the average interaction for each interaction across all datasets.

Arguments: 
1. Pathway to the network dataset (output of the previous script)
2. Pathway for output data

### CD_pycombo.py
CD Application: Following data preparation, we executed the CD analysis using pycombo and set the parameters as "modularity_resolution=1.4" and "max_communities=46" (line 31) so the final result clusters nodes in 46 communities to resemble diploid (2n) genomes. 

Arguments: 
1. pathway to the average network data (output of the previous step)
2. Pathway for output data

## Visualization
The results obtained through CD analysis can be effectively visualized using the following tools:

### Gephi (https://gephi.org/)
Gephi provides a versatile platform for visualizing and analyzing complex networks, making it an ideal choice for exploring the outcomes of CD in 'genome topology maps'.

### Helios web (https://github.com/filipinascimento/helios-web/)
Helios Web offers a web-based solution for visualizing and interpreting the results of CD in 3D.
