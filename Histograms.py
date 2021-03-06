
# coding: utf-8

# In[7]:


import pandas as pd 
import numpy as np
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
os.chdir("./BI_spring_2018/Project")
df =  pd.read_csv('./34_CpGs_10x.csv', sep="\t")



df = df.drop("Unnamed: 0", axis = 1)
df = df.transpose()
df.columns = df.loc["P"]
df = df.drop("P", axis = 0)


# In[11]:


Range = []
Ns = []
with tqdm(total=2611773) as pbar:
    
        
    for position in df:
        pbar.update(1)
        y1 = df[position][df[position] != "-"].astype(np.float)
        Range += [max(y1)-min(y1)]
        Ns += [len(y1)]


# In[13]:




plt.hist(Range, bins = 20, color = "navy")
plt.xlabel('Max-Min')
plt.ylabel('Frequency')
plt.title('Histogram of methylation range distribution')
plt.show()


plt.hist(Ns, color = "darkblue", bins = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34])
plt.xlabel('Number of samples')
plt.ylabel('Count')
plt.title('Histogram samples number')
plt.axis([0, 35, 0, 310000])
plt.show()

