
# coding: utf-8

# In[390]:


import pandas as pd

id_list = pd.read_csv('./RRBS_TCDD_34_selected_subjects.csv', sep=";")
id_list = id_list.drop("id", axis = 1)
id_list.columns = ("ID", "TCDD", "TCDD_quartiles")
id_list = id_list.set_index(id_list["ID"])


# In[170]:
df =  pd.read_csv('./34_CpGs_10x.csv', sep="\t")


# In[171]:


df = df.drop("Unnamed: 0", axis = 1)

df = df.transpose()
df = df.drop("ID", axis = 1)
df.columns = df.loc["P"]
df = df.drop("P")




# In[346]:


import time
import scipy.stats
import numpy as np



x = id_list.transpose()
x = x.drop("ID", axis = 1)
x = x.drop("P")


# In[424]:


#регрессии
slopes = []
intercepts = []
rs = []
ps = []
st_errs = []
means = []
st_ds = []
mins = []
sts_1 = []
medians = []
sts_3= []
maxs = []
i = 0
for position in df:
    y = df[position]
    current_df = pd.concat([x, y], axis = 1)
    current_df = current_df[current_df[position] != "-"]
    y = np.asarray(current_df[position])
    X = np.asarray(current_df["Dioxine"])
    y = y.astype(np.float)
    X = X.astype(np.float)
    if len(y) >= 3:
        slope, intercept, r, p, std_err = scipy.stats.linregress(X, y)
        slopes += [slope]
        intercepts += [intercept]
        rs += [r]
        ps += [p]
        st_errs += [std_err]
        means += [np.mean(y)]
        st_ds += [np.std(y)]
        mins += [min(y)]
        sts_1 += [np.percentile(y, 25)]
        medians += [np.median(y)]
        sts_3 += [np.percentile(y, 75)]
        maxs += [max(y)]
    else:
        slopes += ["-"]
        intercepts += ["-"]
        rs += ["-"]
        ps += ["-"]
        st_errs += ["-"]
        means += ["-"]
        st_ds += ["-"]
        mins += [min(y)]
        sts_1 += ["-"]
        medians += ["-"]
        sts_3 += ["-"]
        maxs += [max(y)]
df.loc["TCDD_B"] = slopes
df.loc["TCDD_A"] = intercepts
df.loc["TCDD_P"] = ps
df.loc["TCDD_R2"] = rs
df.loc["st_err"] = st_errs

df.loc["Mean"] = means
df.loc["St_d"] = st_ds
df.loc["Min"] = mins
df.loc["1st_q"] = sts_1
df.loc["Median"] = medians
df.loc["3rd_q"] = sts_3
df.loc["Max"] = maxs



df_to_csv = df.transpose()



df_to_csv.to_csv("./34_all_CpGs_10x_with_regression_and_CpGs_info.csv", sep = "\t")



