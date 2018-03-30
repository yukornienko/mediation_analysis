
# coding: utf-8

# In[63]:


import pandas as pd
ID_list = ["P"]
i = 0
#читаем и сохраняем список интересующих нас id
with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/RRBS_TCDD_34_selected_subjects.csv") as IDs: 
    for line in IDs:
        tab = line.strip().split(";")
        if i == 0:
            i+=1
        else:
            ID_list += [tab[1]]
#print( ID_list) 
        


# In[64]:


#читаем файл со всеми CpGs с покрытием >10
all_CpGs =  pd.read_csv('/Users/yukornienko/Downloads/BI_spring_2018/Project/all_cpgs_10x.txt', sep="\t")


# In[65]:


our_CpGs = all_CpGs[ID_list]
#our_CpGs - в изначальном файле - 7 132 072 CpGs (без удаления пустых строк).


# In[70]:


#записываем файл с CpGs, присутствующими в хотя бы одном образце, только для итересующих нас id
our_CpGs.to_csv("/Users/yukornienko/Downloads/BI_spring_2018/Project/our_CpGs.csv", sep = "\t")

n = 0
with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/our_CpGs.csv") as data: 
    with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/34_CpGs_10x.csv", "w") as out: 
        for line in data:
            tab = line.strip().split("\t")
            flag = 0
            for i in range(len(tab)):
                if tab[i] == "-":
                    flag += 1
                else:
                    continue
            if flag <= 33:
                out.write(line)
                n += 1
print(n) #n = 2611774 строк => число CpGs = 2 611 773 (т.к. первая строка - названия колонок) 


# In[71]:


#создаем и записываем файл с CpGs, присутствующих в каждом образце среди итересующих нас id
n = 0
our_CpGs_presented_everywhere = pd.DataFrame()
with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/our_CpGs.csv") as data: 
    with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/34_CpGs_10x_present_everywhere.csv", "w") as out: 
        for line in data:
            tab = line.strip().split("\t")
            flag = 0
            for i in range(len(tab)):
                if tab[i] == "-":
                    flag = 1
                else:
                    continue
            if flag == 0:
                out.write(line)
                n += 1
print(n) #n = 29162 строк => число CpGs = 29 161 (т.к. первая строка - названия колонок) 


            


# In[74]:


#фильтруем исходный файл, создавая новый файл без пустых строк
n = 0
with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/all_cpgs_10x.txt") as data: 
    with open ("/Users/yukornienko/Downloads/BI_spring_2018/Project/all_cpgs_10x_filtered.csv", "w") as out: 
        for line in data:
            tab = line.strip().split("\t")
            flag = 0
            for i in range(len(tab)):
                if tab[i] == "-":
                    flag += 1
                else:
                    continue
            if flag <= 39:
                out.write(line)
                n += 1
print(n) # n = 2647798 => в исходном файле для 40 образцов, лишь для 2 647 787 CpGs покрытие было >10X (строки непустые). 

