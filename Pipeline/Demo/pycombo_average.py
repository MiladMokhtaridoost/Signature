# -*- coding: utf-8 -*-
"""
"""
import networkx as nx
import pandas as pd
import pycombo

cells =  [
  "Astrocyte_Spine",
  "H9hESC_day00_Zhang"]

### set the directory where your data is located

for i in range(1, len(cells)+1):
  cell = cells[i-1]
  print(f'{i} = {cell}')
### load the dataset
  data = pd.read_csv(f'{data_folder}/average_1MB_female_network.txt', usecols=['ID_chrA', 'ID_chrB', 'average'], delimiter=' ')

### change the current working directory to the data folder
#os.chdir(data_folder)

#for i in range(1, len(cells)+1):
#for i in range(1, 3):
    
#    cell = cells[i-1]
#    print(f'{i} = {cell}')
### load the dataset
#cell = "average_1MB"
#data = pd.read_csv(f'{data_folder}/average_netwrok_filtered_by_signature.csv', usecols=['ID_chrA', 'ID_chrB', 'freq'], delimiter=' ')
#data = pd.read_csv(f'{data_folder}/average_1MB_male_network.txt', usecols=['ID_chrA', 'ID_chrB', 'average'], delimiter=' ')

data = pd.read_csv(
    f'{data_folder}/average_1MB_male_network.txt',
    usecols=['ID_chrA', 'ID_chrB', 'average'],
    delimiter='\t'
)

print(data.head)
#data = data.head(7000)

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
    
###### Each community in a row      
##with open(f'{cell}_comms.txt', 'w') as file:
    # iterate over each element in the list
##    for com, nodes in comms.items():
        # write the element to the file followed by a newline character
##        file.write(",".join(nodes) + '\n')



##### each node in a row
with open(f'{data_folder}/average_netwrok_male_comms_46.txt', 'w') as file:
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
        
        
        
