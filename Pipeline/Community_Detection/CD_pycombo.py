# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 15:37:26 2023

@author: Milad Mokhtaridoost
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 19:06:11 2023

@author: Milad Mokhtaridoost
"""
import networkx as nx
import pandas as pd
#import os
import pycombo
import sys

data_folder = sys.argv[1]
print(f'data folder = {data_folder}')
result_folder = sys.argv[2]
print(f'result folder = {result_folder}')

### load the dataset
data = pd.read_csv(f'{data_folder}/average_1MB_network.txt', usecols=['ID_chrA', 'ID_chrB', 'freq'], delimiter=' ')

### shifting all weights to positive (in case there is any negative value)
min_weight = data['freq'].min()
data['freq'] = data['freq'] + abs(min_weight)


### constructing network
graph = nx.from_pandas_edgelist(data, source='ID_chrA', target='ID_chrB', edge_attr='freq')

print(nx.number_of_nodes(graph))
print(nx.number_of_edges(graph))

### apply pycombo (46 communities)
communities = pycombo.execute(graph, 'freq', modularity_resolution=1.4, max_communities=46)

comms = {}
for node, c in communities[0].items():
    if c not in comms:
        comms[c] = [node]
    else:
            
        comms[c].append(node)
print(len(comms))
    
##### each node in a row
with open(f'{result_folder}/final_comms.txt', 'w') as file:
    # create a list to store the rows of the dataframe
    rows = []
    # initialize the community ID counter to 1
    com_id = 1
    # iterate over each element in the list
    for com, nodes in comms.items():
    # iterate over each node in the community
        for node in nodes:
            # append a row to the list with the node ID and community ID
            rows.append([node, com_id])
        # increment the community ID counter for the next community
        com_id += 1
    # create a dataframe from the rows list
    df = pd.DataFrame(rows, columns=['ID name', 'Community ID'])
    # write the dataframe to the file as a CSV
    df.to_csv(file, index=False)
