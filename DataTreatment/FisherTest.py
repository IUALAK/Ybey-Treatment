import scipy.stats as stats
import pandas as pd

data = pd.read_csv('./info_table.csv')
data = data.groupby(['ybey_group', 's11_group']).count()
data = [[data.iloc[0,0],data.iloc[1,0]],[data.iloc[2,0],data.iloc[3,0]]]
print(data)
print(stats.fisher_exact(data))