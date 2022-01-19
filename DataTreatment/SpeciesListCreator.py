from Bio import Entrez
import pandas as pd
import re


def species_treater(answer):
    answer = answer['DocumentSummarySet']['DocumentSummary']
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

# pick the needed
requests = [
'"Eukaryota"[Organism] AND (latest[filter] AND ("complete genome"[filter] OR "chromosome level"[filter]) AND (all[filter] NOT anomalous[filter] AND all[filter] NOT partial[filter]))',
# '"Saccharibacteria"[Organism] AND (latest[filter] AND "complete genome"[filter] AND all[filter] NOT anomalous[filter])'            
]

species_list = [] # list for assemblies
for request in requests:
    handle = Entrez.esearch(db = 'assembly', term = request, retmax=8000)
    response = Entrez.read(handle)
    for id in response['IdList']:
        try:
            handle = Entrez.esummary(db = 'assembly', id = id, retmode = 'xml')
        except:
            continue
        answer = Entrez.read(handle, validate = False)
        species_treater(answer)
species_series = pd.Series(species_list)
species_series.to_csv('species_list.csv', index = False)