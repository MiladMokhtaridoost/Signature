# Estimating spatial genome topology with Community Detection

We used Community Detection (CD) to explore spatial genome topology. CD generates clusters of cis and trans interaction frequencies, determines multi-way interactions, and visualizes the CD-derived results in genome topology maps.

#### Requierements
- python/3.8.0
  
  - Pycombo (https://pypi.org/project/pycombo/) 

- R/4.2.1


#### Input Data

To perform CD, we utilized the average interaction frequency (from cooler) of our compendium of 62 individual Hi-C datasets (as detailed in the paper) for each interaction.

#### Utilizing Combo Algorithm

We implemented the Combo algorithm, which is available through the pycombo Python package (https://pycombo.readthedocs.io/en/latest/), as a non-overlapping CD algorithm. This approach comprises two fundamental steps:

## Scripts

### data_processing.R

Data Preparation: In this phase, we prepared the data by converting the interaction frequencies (from cooler) for each interaction across all datasets into a suitable format for subsequent CD analysis.

Arguments: 
1. Index of for loopr (for parallel computing) which should be pass from shell script ($1) to the R script
2. Pathway to the dataset (folder that includes all datasets name with Hi-C data in them)
3. Pathway for output data
4. Resolution of Hi-C data (1000000 in this study)

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
