# Get Data
from Bio import Entrez
import pandas as pd
import re

def search_ybey(species):
    research_table = ybey_table.loc[ybey_table['Organism'].str.contains(species)]
    return 1 if len(research_table) != 0 else 0 

def us11_analysis(sequence):
    sequence = sequence[-30:]
    # End of the new segment
    scan = re.findall('(?=(d.{8}))', sequence)
    if scan != []:
        # print(scan)
        if scan[-1][7:] in ['ng','dg']: #and len(scan) == 1: 
            return 1 # ng and dg groups
        else:
            return 0 # other groups
    else:
        return 'NotRepresentative'
# Preparing list
 
ybey_table = pd.read_csv('Ybey.tab', sep = '\t')
s11_table = pd.read_csv('s11.csv')
taxon_table = pd.read_csv('taxonid.csv')
info_dict = {'species': [],'taxon': [], 'ybey_group': [], 's11_sequence': [], 's11_group': [], 'id': []}
species_list = pd.read_csv('./species_list.csv', index_col=False, squeeze=True) # list for assemblies

# Getting main information
for species in species_list:
    local_s11_table = s11_table.loc[s11_table['description'].str.contains(f'OS={species}')].reset_index(drop = True)
    if len(local_s11_table) == 0:
        continue
    # print(local_s11_table.iloc[0])
    for i in range(len(local_s11_table)):
        result_protein = local_s11_table.iloc[i]['sequence']
        if len(result_protein) >= 100:
            protein_id = local_s11_table.iloc[i]['id']
            s11_group = us11_analysis(result_protein)
            break
    # Micromonas commoda to verify
    try:
        taxon = taxon_table[taxon_table['ScientificName'].str.contains(f'{species}')]['Lineage'].reset_index(drop = True)[0]
    except:
        print('taxon error:', species )
        continue
    ybey_presence = search_ybey(species)
    if (s11_group == 0 or s11_group == 1) and len(taxon.split('; ')) >= 5:
        info_dict['species'].append(species)
        info_dict['taxon'].append(taxon)
        info_dict['ybey_group'].append(ybey_presence)
        info_dict['s11_sequence'].append(result_protein)
        info_dict['s11_group'].append(s11_group)
        info_dict['id'].append(protein_id)
        # Part that should be changed for pipline 
info_table = pd.DataFrame.from_dict(info_dict)
info_table.to_csv('info_table.csv', index = False)
