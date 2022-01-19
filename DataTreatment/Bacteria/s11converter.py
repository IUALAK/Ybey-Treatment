import pandas as pd
import numpy as np
from Bio import SeqIO

fast = SeqIO.parse('s11.fasta', 'fasta')
fast_dict = {'description': [], 'sequence': [], 'id': []}
for i in fast:
    fast_dict['description'].append(i.description)
    fast_dict['id'].append(i.id)
    fast_dict['sequence'].append(str(i.seq).lower())

fast_dict = pd.DataFrame.from_dict(fast_dict)
fast_dict.to_csv('s11.csv', index = False)
