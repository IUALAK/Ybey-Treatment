# First version for S11 based version

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

info_table = pd.read_csv('../DataTreatment/info_table.csv')
# getting visualization by ybey group
fig, axes = plt.subplots(1, 2)
neg_ybey = info_table.loc[info_table['ybey_group'] == 0]
pos_ybey = info_table.loc[info_table['ybey_group'] == 1]
pos_ybey.groupby(['s11_group']).size().plot.bar(subplots = True, ax = axes[0], sharex = False, sharey = False, title = 'uS11m distrubution in taxons by Ybey group')
axes[0].set_title('Ybey +')
neg_ybey.groupby(['s11_group']).size().plot.bar(subplots = True, ax = axes[1], sharex = False, sharey = False)
axes[1].set_title('Ybey -')
plt.savefig('./Plots/ybey_plot.png')

# getting hex plots
for domain in ['Eukaryota','Bacteria']:
    ddf = info_table.loc[info_table['taxon'].str.contains(domain)]
    if len(ddf) == 0:
        continue
    ddf = ddf.drop(columns = 'id')
    plot = ddf.plot.hexbin(x = 's11_group', y = 'ybey_group', gridsize = (2, 1), title = f'distribution for {domain}', sharex = False, xticks = (0, 1), yticks = (0, 1) )
    corr = ddf.corr()
    print(corr)
    plot.get_figure().savefig(f'./Plots/{domain}_plot.png')
    step = 0
    for genus in set(ddf['taxon']):
        # family = genus.split('; ')[:-2]
        genus_title = genus.split('; ')[-1]
        gdf = ddf.loc[ddf['taxon'].str.contains(genus)]
        if len(gdf) < 20:
            continue
        plot = gdf.plot.hexbin(x = 's11_group', y = 'ybey_group', gridsize = (1,2), title = f'distribution for {genus_title}', sharex = False, xticks = (0, 1), yticks = (0, 1),extent = (0.7,1,0,1) )
        plot.get_figure().savefig(f'./Plots/{genus_title}_plot.png')
        step += 1
        if step == 10:
            step = 0
            break


