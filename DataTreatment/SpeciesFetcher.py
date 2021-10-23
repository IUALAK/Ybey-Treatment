from Bio import Entrez
import pandas as pd
import re

def species_treater(answer):
    answer = answer['DocumentSummarySet']['DocumentSummary']
    # for i in answer:
    #    answer =  i['Organism'] #['SpeciesName']
    try:
        answer = answer[0]['SpeciesName']
    except:
        return 0
    if re.search(r'(\bbacterium\b)$|(\bsp\.)$', answer) != None:
        return 0
    answer = re.sub('\[|\]', '', answer)
    if answer not in species_list:
        species_list.append(answer)
        # print(answer)
    return 0
 

Entrez.email = 'alakhrasyunes@gmail.com'

requests = [
            'Bacteria[filter] AND latest[filter] AND "complete genome"[filter] NOT anomalous[filter]',
            # '(Animals[filter] OR Plants[filter] OR Fungi[filter] OR Protists[filter]) AND latest[filter] AND "complete genome"[filter] NOT anomalous[filter]',
# '"Saccharibacteria"[Organism] AND (latest[filter] AND "complete genome"[filter] AND all[filter] NOT anomalous[filter])'            
]
ybey_table = pd.read_csv('Ybey.tab', sep = '\t')
s11_table = pd.read_csv('s11.csv')
taxon_table = pd.read_csv('taxid.tab', sep = '\t')
info_dict = {'species': [],'taxon': [], 'ybey_group': [], 's11_sequence': [], 's11_group': [], 'id': []}
species_list = [] # list for assemblies
for request in requests:
    handle = Entrez.esearch(db = 'assembly', term = request, retmax=40)
    response = Entrez.read(handle)
    for id in response['IdList']:
        handle = Entrez.esummary(db = 'assembly', id = id, retmode = 'xml')
        answer = Entrez.read(handle, validate = False)
        species_treater(answer)
species_series = pd.Series(species_list)
species_series.to_csv('species_list.csv', index = False)
print('species list is ready for using')