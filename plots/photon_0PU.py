import os
from glob import glob
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib


# ## Configuration

# ### Input files and FE algorithms to be studied
# 0PU photon data are loaded

# In[5]:


algo_files = {}
directory = 'Documents/FPGAs/'
fes = ['Threshold']

for fe in fes:
    algo_files[fe] = [ os.path.join(os.environ['HOME'], directory, 'out.hdf5') ]


# ## Loading and preprocessing dataframes

# In[6]:


algos_dfs = {}
for fe,files in algo_files.items():
    name = fe
    dfs = []
    for file in files:
        with pd.HDFStore(file, mode='r') as store:
            dfs.append(store[name])
    algos_dfs[fe] = pd.concat(dfs)
algo_names = sorted(algos_dfs.keys())


# In[7]:


#Cleaning
df = algos_dfs['Threshold']
df = df[ df['genpart_exeta']>0 ]
df = df[ df['cl3d_eta']>0 ]


# In[8]:


df['enres'] = ( df['genpart_energy']-df['cl3d_energy'] ) / df['genpart_energy']


# In[9]:


nansel = pd.isna(df['enres']) 
nandf = df[nansel]
nandf['enres'] = 1.1
df = df[~nansel]
df = pd.concat([df,nandf], sort=False)
counts, edges = np.histogram(df['enres'].array, bins=500)
nandf


# In[10]:


matplotlib.rcParams.update({'font.size': 11})
def plot_histo_quick(df1, df2, var, nbins=600):
    plt.figure(figsize=(12,4))
    plt.subplot(211)
    plt.hist(df1[var], bins=nbins)
    plt.grid()
    plt.xlabel(var)
    plt.ylabel('Counts')
    plt.subplot(221)
    plt.hist(df2[var], bins=nbins)
    plt.grid()
    plt.xlabel(var)
    plt.ylabel('Counts')
    plt.show()


# In[11]:


enrescut = 0.3
duplicated = df[ df['enres'] > enrescut ]
duplicated.columns


# In[ ]:


for col in duplicated.columns:
    if 'cl3d' in col:
        plot_histo_quick(df, duplicated, var=col, nbins=200)


# In[ ]:




