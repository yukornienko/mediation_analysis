
# coding: utf-8

# In[24]:


import pandas as pd 
import numpy as np
import statsmodels.formula.api as sm
from tqdm import tqdm
import os
os.chdir("/Users/yukornienko/Downloads/BI_spring_2018/Project")
df =  pd.read_csv('./34_CpGs_10x.csv', sep="\t")
print("Hello")


df = df.drop("Unnamed: 0", axis = 1)
df = df.transpose()
df.columns = df.loc["P"]
positions = [df.loc["P"]]
df = df.drop("P", axis = 0)

#df = df.transpose()[1958850:]
#df = df.transpose()[0:1000]
df = df.convert_objects(convert_numeric=True)



# In[25]:


id_list = pd.read_csv('./RRBS_TCDD_34_selected_subjects.csv', sep=";")
id_list = id_list.drop("id", axis = 1)
id_list.columns = ("ID", "TCDD", "TCDD_quartiles")
id_list = id_list.set_index(id_list["ID"])
#print(id_list)

smoking = pd.read_csv('./Smoke_dummies.csv', sep=",")
smoking = smoking.set_index(id_list["ID"])
#print(smoking)


x1 = pd.DataFrame(data = {"TCDD": id_list["TCDD"],"smoking6": smoking["smoking_last6months"]})
x1["TCDD"] = x1["TCDD"].astype(np.float)
x1["smoking6"] = x1["smoking6"].astype('category')

x1 = x1.set_index(id_list["ID"])
x2 = pd.get_dummies(x1)
x2 = x2.set_index(id_list["ID"])
x2["smoking6_0"] = x2["smoking6_0"].astype('category')
x2["smoking6_1"] = x2["smoking6_1"].astype('category')
x2["smoking6_2"] = x2["smoking6_2"].astype('category')
x2["smoking6_3"] = x2["smoking6_3"].astype('category')
x2["smoking6_4"] = x2["smoking6_4"].astype('category')
print(x2)
df = df.set_index(id_list["ID"])
df


# In[10]:


#smoke6 and TCDD conc
intercepts = []
B_S_T0s = []
B_S_T1s = []
B_S_T2s = []
B_S_T3s = []
B_S_T4s = []
B_TCDDs = [] 
rs = []
p_intercepts = []
p_S_T0s = []
p_S_T1s = []
p_S_T2s = []
p_S_T3s = []
p_S_T4s = []
p_TCDDs = []
p_summs = []
means = []
st_ds = []
mins = []
sts_1 = []
medians = []
sts_3= []
maxs = []
Ns = []

#df = df
print("I am starting")

with tqdm(total=2611773) as pbar:
    
        
    for position in df:
        pbar.update(1)
        y1 = df[position]
        current_df = pd.concat([x2, y1], axis = 1)
        current_df.columns = ["TCDD", "smoke0", "smoke1", "smoke2", "smoke3", "smoke4", "position"]
        y2 = y1[np.logical_not(np.isnan(y1))]
        if (max(y2) - min(y2)) >= 20 and len(y2) >= 10:
            
            model=sm.ols(formula = 'position ~ TCDD + smoke0 + smoke1 + smoke2 + smoke3 +smoke4', data = current_df, missing="drop").fit()
            
            rs += [model.rsquared]
            p_intercepts += [model.pvalues[0]]
            intercepts += [model.params[0]]
            p_summs += [model.f_pvalue]
            B_TCDDs += [model.params[6]]
            p_TCDDs += [model.pvalues[6]]
            p_S_T0s += [model.pvalues[1]]
            p_S_T1s += [model.pvalues[2]]
            p_S_T2s += [model.pvalues[3]]
            p_S_T3s += [model.pvalues[4]]                
            p_S_T4s += [model.pvalues[5]]
            B_S_T0s += [model.params[1]]
            B_S_T1s += [model.params[2]]
            B_S_T2s += [model.params[3]]
            B_S_T3s += [model.params[4]]
            B_S_T4s += [model.params[5]]
        else:
            rs += ["-"]
            p_intercepts += ["-"]
            intercepts += ["-"]
            p_summs += ["-"]
            B_TCDDs += ["-"]
            p_TCDDs += ["-"]
            p_S_T0s += ["-"]
            p_S_T1s += ["-"]
            p_S_T2s += ["-"]
            p_S_T3s += ["-"]                
            p_S_T4s += ["-"]
            B_S_T0s += ["-"]
            B_S_T1s += ["-"]
            B_S_T2s += ["-"]
            B_S_T3s += ["-"]
            B_S_T4s += ["-"]

        means += [np.mean(y2)]
        st_ds += [np.std(y2)]
        mins += [min(y2)]
        sts_1 += [np.percentile(y2, 25)]
        medians += [np.median(y2)]
        sts_3 += [np.percentile(y2, 75)]
        maxs += [max(y2)]

        Ns += [len(y2)]




# In[22]:


new_df = pd.DataFrame(data = {"A": intercepts, "B_Smoke_T0": B_S_T0s, "B_Smoke_T1": B_S_T1s, "B_smoke_T2": B_S_T2s, "B_smoke_T3": B_S_T3s, "B_smoke_T4": B_S_T4s, "B_TCDD": B_TCDDs, "R2": rs, "P_A": p_intercepts, "P_smoke_T0": p_S_T0s, "P_smoke_T1": p_S_T1s, "p_smoke_T2": p_S_T2s, "P_smoke_T3": p_S_T3s, "P_smoke_T4":  p_S_T4s,"P_TCDD": p_TCDDs, "P_summ": p_summs, "Mean": means, "St_d": st_ds, "Min": mins, "1st_q": sts_1, "Median": medians, "3rd_q": sts_3, "Max": maxs, "N": Ns})

new_df = new_df.transpose()
new_df.columns = df1.loc["P"]



# In[28]:


new_df = new_df.transpose()
new_df = new_df[new_df["R2"] != "-"]
new_df.to_csv("./34_10x_CpG~TCDD_and_Smoke6_cat.csv", sep = "\t")


# In[29]:


#smoking6
intercepts = []
B_S_T0s = []
B_S_T1s = []
B_S_T2s = []
B_S_T3s = []
B_S_T4s = []
B_TCDDs = [] 
rs = []
p_intercepts = []
p_S_T0s = []
p_S_T1s = []
p_S_T2s = []
p_S_T3s = []
p_S_T4s = []
p_TCDDs = []
p_summs = []
means = []
st_ds = []
mins = []
sts_1 = []
medians = []
sts_3= []
maxs = []
Ns = []

#df = df
print("I am starting")

with tqdm(total=2611773) as pbar:
    
        
    for position in df:
        pbar.update(1)
        y1 = df[position]
        current_df = pd.concat([x2, y1], axis = 1)
        current_df.columns = ["TCDD", "smoke0", "smoke1", "smoke2", "smoke3", "smoke4", "position"]
        y2 = y1[np.logical_not(np.isnan(y1))]
        if (max(y2) - min(y2)) >= 20 and len(y2) >= 10:
            
            model=sm.ols(formula = 'position ~ smoke0 + smoke1 + smoke2 + smoke3 +smoke4', data = current_df, missing="drop").fit()
            
            rs += [model.rsquared]
            p_intercepts += [model.pvalues[0]]
            intercepts += [model.params[0]]
            p_summs += [model.f_pvalue]
          #  B_TCDDs += [model.params[1]]
          #  p_TCDDs += [model.pvalues[1]]
            p_S_T0s += [model.pvalues[1]]
            p_S_T1s += [model.pvalues[2]]
            p_S_T2s += [model.pvalues[3]]
            p_S_T3s += [model.pvalues[4]]                
            p_S_T4s += [model.pvalues[5]]
            B_S_T0s += [model.params[1]]
            B_S_T1s += [model.params[2]]
            B_S_T2s += [model.params[3]]
            B_S_T3s += [model.params[4]]
            B_S_T4s += [model.params[5]]
        else:
            rs += ["-"]
            p_intercepts += ["-"]
            intercepts += ["-"]
            p_summs += ["-"]
           # B_TCDDs += ["-"]
           # p_TCDDs += ["-"]
            p_S_T0s += ["-"]
            p_S_T1s += ["-"]
            p_S_T2s += ["-"]
            p_S_T3s += ["-"]                
            p_S_T4s += ["-"]
            B_S_T0s += ["-"]
            B_S_T1s += ["-"]
            B_S_T2s += ["-"]
            B_S_T3s += ["-"]
            B_S_T4s += ["-"]

        means += [np.mean(y2)]
        st_ds += [np.std(y2)]
        mins += [min(y2)]
        sts_1 += [np.percentile(y2, 25)]
        medians += [np.median(y2)]
        sts_3 += [np.percentile(y2, 75)]
        maxs += [max(y2)]

        Ns += [len(y2)]




# In[30]:


new_df = pd.DataFrame(data = {"A": intercepts, "B_Smoke_T0": B_S_T0s, "B_Smoke_T1": B_S_T1s, "B_Smoke_T2": B_S_T2s, "B_Smoke_T3": B_S_T3s, "B_Smoke_T4": B_S_T4s, "R2": rs, "P_A": p_intercepts, "P_Smoke_T0": p_S_T0s, "P_Smoke_T1": p_S_T1s, "P_Smoke_T2": p_S_T2s, "P_Smoke_T3": p_S_T3s, "P_Smoke_T4":  p_S_T4s, "P_summ": p_summs, "Mean": means, "St_d": st_ds, "Min": mins, "1st_q": sts_1, "Median": medians, "3rd_q": sts_3, "Max": maxs, "N": Ns})

new_df = new_df.transpose()
new_df.columns = positions
new_df = new_df.transpose()
new_df = new_df[new_df["R2"] != "-"]
new_df.to_csv("./34_10x_CpG~Smoke6_cat.csv", sep = "\t")



# In[35]:


id_list = pd.read_csv('./RRBS_TCDD_34_selected_subjects.csv', sep=";")
id_list = id_list.drop("id", axis = 1)
id_list.columns = ("ID", "TCDD", "TCDD_quartiles")
id_list = id_list.set_index(id_list["ID"])
#print(id_list)

smoking = pd.read_csv('./Smoke_dummies.csv', sep=",")
smoking = smoking.set_index(id_list["ID"])
#print(smoking)


x1 = pd.DataFrame(data = {"TCDD": id_list["TCDD_quartiles"],"smoking3": smoking["smoking_last6months_reduced"]})
x1["TCDD"] = x1["TCDD"].astype('category')
x1["smoking3"] = x1["smoking3"].astype('category')

x1 = x1.set_index(id_list["ID"])
x2 = pd.get_dummies(x1)
print(x1)
x2 = x2.set_index(id_list["ID"])
print(x2)
x2["TCDD_1"] = x2["TCDD_1"].astype('category')
x2["TCDD_2"] = x2["TCDD_2"].astype('category')
x2["TCDD_3"] = x2["TCDD_3"].astype('category')
x2["smoking3_0"] = x2["smoking3_0"].astype('category')
x2["smoking3_1"] = x2["smoking3_1"].astype('category')
x2["smoking3_2"] = x2["smoking3_2"].astype('category')
#x2["smoking6_3"] = x2["smoking6_3"].astype('category')
#x2["smoking6_4"] = x2["smoking6_4"].astype('category')
print(x2)



# In[37]:


df


# In[39]:


###smoking3
intercepts = []
B_S_T0s = []
B_S_T1s = []
B_S_T2s = []
B_S_T3s = []
B_S_T4s = []
B_TCDDs = [] 
rs = []
p_intercepts = []
p_S_T0s = []
p_S_T1s = []
p_S_T2s = []
p_S_T3s = []
p_S_T4s = []
p_TCDDs = []
p_summs = []
means = []
st_ds = []
mins = []
sts_1 = []
medians = []
sts_3= []
maxs = []
Ns = []

#df = df
print("I am starting")
i = 0
with tqdm(total=2611773) as pbar:
#with tqdm(total=100) as pbar: 
        
    for position in df:
        pbar.update(1)
        y1 = df[position]
        current_df = pd.concat([x2, y1], axis = 1)
        current_df.columns = ["TCDD1", "TCDD2", "TCDD3", "smoke0", "smoke1", "smoke2", "position"]
        y2 = y1[np.logical_not(np.isnan(y1))]
        if (max(y2) - min(y2)) >= 20 and len(y2) >= 10:
            
            model=sm.ols(formula = 'position ~ smoke0 + smoke1 + smoke2', data = current_df, missing="drop").fit()
            
            rs += [model.rsquared]
            p_intercepts += [model.pvalues[0]]
            intercepts += [model.params[0]]
            p_summs += [model.f_pvalue]
          #  B_TCDDs += [model.params[1]]
          #  p_TCDDs += [model.pvalues[1]]
            p_S_T0s += [model.pvalues[1]]
            p_S_T1s += [model.pvalues[2]]
            p_S_T2s += [model.pvalues[3]]
            #p_S_T3s += [model.pvalues[4]]                
            #p_S_T4s += [model.pvalues[5]]
            B_S_T0s += [model.params[1]]
            B_S_T1s += [model.params[2]]
            B_S_T2s += [model.params[3]]
            #B_S_T3s += [model.params[4]]
            #B_S_T4s += [model.params[5]]
           # print(model.summary())
        else:
            rs += ["-"]
            p_intercepts += ["-"]
            intercepts += ["-"]
            p_summs += ["-"]
           # B_TCDDs += ["-"]
           # p_TCDDs += ["-"]
            p_S_T0s += ["-"]
            p_S_T1s += ["-"]
            p_S_T2s += ["-"]
           # p_S_T3s += ["-"]                
           # p_S_T4s += ["-"]
            B_S_T0s += ["-"]
            B_S_T1s += ["-"]
            B_S_T2s += ["-"]
            #B_S_T3s += ["-"]
            #B_S_T4s += ["-"]

        means += [np.mean(y2)]
        st_ds += [np.std(y2)]
        mins += [min(y2)]
        sts_1 += [np.percentile(y2, 25)]
        medians += [np.median(y2)]
        sts_3 += [np.percentile(y2, 75)]
        maxs += [max(y2)]

        Ns += [len(y2)]


new_df = pd.DataFrame(data = {"A": intercepts, "B_Smoke_T0": B_S_T0s, "B_Smoke_T1": B_S_T1s, "B_Smoke_T2": B_S_T2s, "R2": rs, "P_A": p_intercepts, "P_Smoke_T0": p_S_T0s, "P_Smoke_T1": p_S_T1s, "P_Smoke_T2": p_S_T2s, "P_summ": p_summs, "Mean": means, "St_d": st_ds, "Min": mins, "1st_q": sts_1, "Median": medians, "3rd_q": sts_3, "Max": maxs, "N": Ns})

new_df = new_df.transpose()
new_df.columns = positions
new_df = new_df.transpose()
new_df = new_df[new_df["R2"] != "-"]
new_df



# In[40]:


new_df.to_csv("./34_10x_CpG~Smoke3_cat.csv", sep = "\t")


# In[41]:


#terciles dioxines
intercepts = []
B_S_T0s = []
B_S_T1s = []
B_S_T2s = []
B_S_T3s = []
B_S_T4s = []
B_TCDDs = [] 
rs = []
p_intercepts = []
p_S_T0s = []
p_S_T1s = []
p_S_T2s = []
p_S_T3s = []
p_S_T4s = []
p_TCDDs = []
p_summs = []
means = []
st_ds = []
mins = []
sts_1 = []
medians = []
sts_3= []
maxs = []
Ns = []
p_TCDD0_s = []
p_TCDD1_s = []
p_TCDD2_s = []
B_TCDD0_s = []
B_TCDD1_s = []
B_TCDD2_s = []

#df = df
print("I am starting")

with tqdm(total=2611773) as pbar:
    
        
    for position in df:
        pbar.update(1)
        y1 = df[position]
        current_df = pd.concat([x2, y1], axis = 1)
        current_df.columns = ["TCDD0", "TCDD1", "TCDD2", "smoke0", "smoke1", "smoke2", "position"]
        y2 = y1[np.logical_not(np.isnan(y1))]
        if (max(y2) - min(y2)) >= 20 and len(y2) >= 10:
            
            model=sm.ols(formula = 'position ~ TCDD0 + TCDD1 + TCDD2', data = current_df, missing="drop").fit()
            
            rs += [model.rsquared]
            p_intercepts += [model.pvalues[0]]
            intercepts += [model.params[0]]
            p_summs += [model.f_pvalue]
          #  B_TCDDs += [model.params[1]]
          #  p_TCDDs += [model.pvalues[1]]
            p_TCDD0_s += [model.pvalues[1]]
            p_TCDD1_s += [model.pvalues[2]]
            p_TCDD2_s += [model.pvalues[3]]
           # p_S_T1s += [model.pvalues[2]]
           # p_S_T2s += [model.pvalues[3]]
            #p_S_T3s += [model.pvalues[4]]                
            #p_S_T4s += [model.pvalues[5]]
            B_TCDD0_s += [model.params[1]]
            B_TCDD1_s += [model.params[2]]
            B_TCDD2_s += [model.params[3]]
            #B_S_T3s += [model.params[4]]
            #B_S_T4s += [model.params[5]]
        else:
            rs += ["-"]
            p_intercepts += ["-"]
            intercepts += ["-"]
            p_summs += ["-"]
           # B_TCDDs += ["-"]
           # p_TCDDs += ["-"]
            p_TCDD0_s += ["-"]
            p_TCDD1_s += ["-"]
            p_TCDD2_s += ["-"]
           # p_S_T3s += ["-"]                
           # p_S_T4s += ["-"]
            B_TCDD0_s += ["-"]
            B_TCDD1_s += ["-"]
            B_TCDD2_s += ["-"]
            #B_S_T3s += ["-"]
            #B_S_T4s += ["-"]

        means += [np.mean(y2)]
        st_ds += [np.std(y2)]
        mins += [min(y2)]
        sts_1 += [np.percentile(y2, 25)]
        medians += [np.median(y2)]
        sts_3 += [np.percentile(y2, 75)]
        maxs += [max(y2)]

        Ns += [len(y2)]


new_df = pd.DataFrame(data = {"A": intercepts, "B_TCDD_T0": B_TCDD0_s, "B_TCDD_T1": B_TCDD1_s, "B_TCDD_T2": B_TCDD2_s, "R2": rs, "P_A": p_intercepts, "P_TCDD_T0": p_TCDD0_s, "P_TCDD_T1": p_TCDD1_s, "P_TCDD_T2": p_TCDD2_s, "P_summ": p_summs, "Mean": means, "St_d": st_ds, "Min": mins, "1st_q": sts_1, "Median": medians, "3rd_q": sts_3, "Max": maxs, "N": Ns})

new_df = new_df.transpose()
new_df.columns = positions
new_df = new_df.transpose()
new_df = new_df[new_df["R2"] != "-"]
new_df.to_csv("./34_10x_CpG~TCDD_terc_cat.csv", sep = "\t")


# In[45]:


#terciles dioxines + smoking3
intercepts = []
B_S_T0s = []
B_S_T1s = []
B_S_T2s = []
B_S_T3s = []
B_S_T4s = []
B_TCDDs = [] 
rs = []
p_intercepts = []
p_S_T0s = []
p_S_T1s = []
p_S_T2s = []
p_S_T3s = []
p_S_T4s = []
p_TCDDs = []
p_summs = []
means = []
st_ds = []
mins = []
sts_1 = []
medians = []
sts_3= []
maxs = []
Ns = []
p_TCDD0_s = []
p_TCDD1_s = []
p_TCDD2_s = []
B_TCDD0_s = []
B_TCDD1_s = []
B_TCDD2_s = []

#df = df
print("I am starting")

with tqdm(total=2611773) as pbar:
    
        
    for position in df:
        pbar.update(1)
        y1 = df[position]
        current_df = pd.concat([x2, y1], axis = 1)
        current_df.columns = ["TCDD1", "TCDD2", "TCDD3", "smoke0", "smoke1", "smoke2", "position"]
        y2 = y1[np.logical_not(np.isnan(y1))]
        if (max(y2) - min(y2)) >= 20 and len(y2) >= 10:
            
            model=sm.ols(formula = 'position ~ TCDD1 + TCDD2 + TCDD3 + smoke0 + smoke1 + smoke2', data = current_df, missing="drop").fit()
            
            rs += [model.rsquared]
            p_intercepts += [model.pvalues[0]]
            intercepts += [model.params[0]]
            p_summs += [model.f_pvalue]
            p_TCDD0_s += [model.pvalues[1]]
            p_TCDD1_s += [model.pvalues[2]]
            p_TCDD2_s += [model.pvalues[3]]
            p_S_T0s += [model.pvalues[4]]
            p_S_T1s += [model.pvalues[5]]
            p_S_T2s += [model.pvalues[6]]                

            B_TCDD0_s += [model.params[1]]
            B_TCDD1_s += [model.params[2]]
            B_TCDD2_s += [model.params[3]]
            B_S_T0s += [model.params[4]]
            B_S_T1s += [model.params[5]]
            B_S_T2s += [model.params[6]]
        else:
            rs += ["-"]
            p_intercepts += ["-"]
            intercepts += ["-"]
            p_summs += ["-"]
      
            p_TCDD0_s += ["-"]
            p_TCDD1_s += ["-"]
            p_TCDD2_s += ["-"]
            p_S_T0s += ["-"]
            p_S_T1s += ["-"]
            p_S_T2s += ["-"]                

            B_TCDD0_s += ["-"]
            B_TCDD1_s += ["-"]
            B_TCDD2_s += ["-"]
            B_S_T0s += ["-"]
            B_S_T1s += ["-"]
            B_S_T2s += ["-"]

        means += [np.mean(y2)]
        st_ds += [np.std(y2)]
        mins += [min(y2)]
        sts_1 += [np.percentile(y2, 25)]
        medians += [np.median(y2)]
        sts_3 += [np.percentile(y2, 75)]
        maxs += [max(y2)]

        Ns += [len(y2)]


new_df = pd.DataFrame(data = {"A": intercepts, "B_Smoke_T0": B_S_T0s, "B_Smoke_T1": B_S_T1s, "B_Smoke_T2": B_S_T2s, "B_TCDD_T1": B_TCDD0_s, "B_TCDD_T2": B_TCDD1_s, "B_TCDD_T3": B_TCDD2_s, "R2": rs, "P_A": p_intercepts, "P_Smoke_T0": p_S_T0s, "P_Smoke_T1": p_S_T1s, "P_Smoke_T2": p_S_T2s, "P_TCDD_T1": p_TCDD0_s, "P_TCDD_T2": p_TCDD1_s, "P_TCDD_T3": p_TCDD2_s, "P_summ": p_summs, "Mean": means, "St_d": st_ds, "Min": mins, "1st_q": sts_1, "Median": medians, "3rd_q": sts_3, "Max": maxs, "N": Ns})

new_df = new_df.transpose()
new_df.columns = positions
new_df = new_df.transpose()
new_df = new_df[new_df["R2"] != "-"]
new_df.to_csv("./34_10x_CpG~TCDD_terc_and_smoke3_cat.csv", sep = "\t")

