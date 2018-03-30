
# coding: utf-8

# In[390]:


import pandas as pd
ID_list = ["P"]
Dioxines = ["Dioxine"]
i = 0
#читаем и сохраняем список интересующих нас id
with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/RRBS_TCDD_34_selected_subjects.csv") as IDs: 
    for line in IDs:
       # print(line)
        tab = line.strip().split(";")
        if i == 0:
            i+=1
        else:
            ID_list += [tab[1]]
            Dioxines += [tab[2]]
#print(ID_list) 
#print (Dioxines)
id_list = pd.DataFrame(data = {"ID": ID_list, 'Dioxine': Dioxines})


# In[391]:


id_list.to_csv("/Users/yukornienko/Downloads/BI_spring_2018/Project/id_list.csv", sep = "\t")


# In[239]:


id_list = id_list.transpose()


# In[240]:


id_list.columns = ID_list
id_list


# In[170]:


our_CpGs =  pd.read_csv('/Users/yukornienko/Downloads/BI_spring_2018/Project/34_CpGs_10x.csv', sep="\t")


# In[171]:


our_CpGs = our_CpGs.drop("Unnamed: 0", axis = 1)


# In[242]:


df = pd.concat([id_list, our_CpGs])


# In[243]:


df = df.transpose()


# In[244]:


df = df.drop("ID", axis = 1)


# In[245]:


df.columns = df.loc["P"]


# In[246]:


df = df.drop("P")


# In[248]:


df.to_csv("/Users/yukornienko/Downloads/BI_spring_2018/Project/our_CpGs_with_dioxines.csv", sep = "\t")


# In[ ]:


###### далее будем строить 2611773 регрессионные модели и отбирать статистически значимые!


# In[346]:


import time
import scipy.stats
import numpy as np


# In[458]:


df


# In[347]:


df_1 = df.drop("Dioxine", axis = 1)


# In[459]:


id_list


# In[351]:


x = id_list.transpose()
x = x.drop("ID", axis = 1)
x = x.drop("P")
print(x)


# In[424]:


#регрессии
slopes = []
intercepts = []
rs = []
ps = []
st_errs = []
for position in df_1:
    y = df_2[position]
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
    else:
        slopes += ["-"]
        intercepts += ["-"]
        rs += ["-"]
        ps += ["-"]
        st_errs += ["-"]
df_1.loc["TCDD_B"] = slopes
df_1.loc["TCDD_A"] = intercepts
df_1.loc["TCDD_P"] = ps
df_1.loc["TCDD_R2"] = rs
df_1.loc["st_err"] = st_errs
df_1
print(time.clock())


# In[428]:


#описательные статистики
df_2 = df_1[:34]
means = []
st_ds = []
mins = []
sts_1 = []
medians = []
sts_3= []
maxs = []
i = 0
for position in df_2:
    i += 1
    y = df_2[position]
    current_df = pd.concat([x, y], axis = 1)
    current_df = current_df[current_df[position] != "-"]
    y = np.asarray(current_df[position])
    y = y.astype(np.float)

    if len(y) >= 3:
        means += [np.mean(y)]
        st_ds += [np.std(y)]
        mins += [min(y)]
        sts_1 += [np.percentile(y, 25)]
        medians += [np.median(y)]
        sts_3 += [np.percentile(y, 75)]
        maxs += [max(y)]
        
    else:
        means += ["-"]
        st_ds += ["-"]
        mins += [min(y)]
        sts_1 += ["-"]
        medians += ["-"]
        sts_3 += ["-"]
        maxs += [max(y)]

df_1.loc["Mean"] = means
df_1.loc["St_d"] = st_ds
df_1.loc["Min"] = mins
df_1.loc["1st_q"] = sts_1
df_1.loc["Median"] = medians
df_1.loc["3rd_q"] = sts_3
df_1.loc["Max"] = maxs
df_1
print(time.clock())


# In[429]:


df_to_csv = df_1.transpose()


# In[466]:


df_to_csv


# In[430]:


df_to_csv.to_csv("/Users/yukornienko/Downloads/BI_spring_2018/Project/34_all_CpGs_10x_with_regression_and_CpGs_info.csv", sep = "\t")


# In[461]:


df_to_csv_2 = df_to_csv[df_to_csv["TCDD_P"] != '-']
df_to_csv_2.to_csv("/Users/yukornienko/Downloads/BI_spring_2018/Project/34_CpGs_10x_with_regression_and_CpGs_info_present_in3.csv", sep = "\t")


# In[467]:


df_to_csv_2


# In[465]:


df_sign = df_to_csv_2[(df_to_csv_2["TCDD_P"]).astype(np.float) < 0.05]
df_sign = df_sign[(abs(df_sign["TCDD_R2"].astype(np.float))) > 0.80]
df_sign = df_sign.transpose()
df_sign.to_csv("/Users/yukornienko/Downloads/BI_spring_2018/Project/34_CpGs_10x_with_regression_and_CpGs_info_sing.csv", sep = "\t")


# In[468]:


df_sign

