# -*- coding: utf-8 -*-
"""
"""
import networkx as nx
import pandas as pd
import pycombo

### set the directory where your data is located
data_folder = '../Demo'

### load the dataset

data = pd.read_csv(
    f'{data_folder}/average_1MB_network.txt',
    usecols=['ID_chrA', 'ID_chrB', 'average'],
    delimiter='\t'
)

print(data.head)

    #if __name__ == "__main__":

min_weight = data['average'].min()
data['average'] = data['average'] + abs(min_weight)

graph = nx.from_pandas_edgelist(data, source='ID_chrA', target='ID_chrB', edge_attr='average')

print(nx.number_of_nodes(graph))
print(nx.number_of_edges(graph))
    # Solve
communities = pycombo.execute(graph, 'average', modularity_resolution=1.2, max_communities=46)

comms = {}
for node, c in communities[0].items():
    if c not in comms:
            comms[c] = [node]
    else:
            
            comms[c].append(node)
print(len(comms))
    
##### each node in a row
with open(f'{data_folder}/average_network_comms.txt', 'w') as file:
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
        
