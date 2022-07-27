#!/usr/bin/env python
# coding: utf-8

# In[1]:


import json
import pandas as pd


# In[2]:


table=pd.read_csv('data_table.tsv', sep='\t')
table


# In[3]:


table["orientation"][0]


# In[4]:


f = open("data_report.jsonl")
data= json.load(f)
#data
#nested dictionaries and lists
exons=data['transcripts'][0]['genomicLocations'][0]['exons']


# In[5]:


exons


# In[17]:


def annotate_gene():
    ann_gene=dict()
    count=1
    s=1
    rna_begin=1
    
    for i in range(len(exons)):
        
        
        if table["orientation"][0]== "+;+":
            #exon position
            begin=int(exons[i]["begin"])
            end=int(exons[i]["end"])



            lenEx = end - begin  + 1

            ann_gene["exon " + str(count)]=[lenEx,(s,s+lenEx-1),(begin,end),(rna_begin, rna_begin+lenEx-1)]
            rna_begin=rna_begin+lenEx
            
            #intron positions
            if (count < len(exons)):
                #intron
                intEnd=int(exons[i+1]["begin"])-1

                intBegin=end+1

                lenInt= intEnd - intBegin + 1

                ann_gene["intron " + str(count)]=[lenInt,(s+lenEx,s+lenEx+lenInt-1),(intBegin,intEnd), "-"]

                s=s+lenEx+lenInt
                count=count+1
                
                
        elif table["orientation"][0]== "-":
            #exon position
            begin=int(exons[i]["end"])
            end=int(exons[i]["begin"])

            lenEx = begin - end  + 1

            ann_gene["exon " + str(count)]=[lenEx,(s,s+lenEx-1),(begin,end),(rna_begin, rna_begin+lenEx-1)]
            rna_begin=rna_begin+lenEx
            

            if (count < len(exons)):
                #intron
                intEnd=int(exons[i+1]["end"])+1
                intBegin=end-1

                lenInt= intBegin - intEnd + 1

                ann_gene["intron " + str(count)]=[lenInt,(s+lenEx,s+lenEx+lenInt-1),(intBegin,intEnd),"-"]

                s=s+lenEx+lenInt
                count=count+1
    return ann_gene


# In[18]:


gene=annotate_gene()


# In[19]:


gene


# In[20]:


df=pd.DataFrame(gene).T
df.columns=["exon/intron len", "position", "pos on chrom", "pos on rna-seq"]


# In[21]:


df


# In[ ]:





# In[ ]:





# In[ ]:





# In[22]:


g = open("gene.fna", "r")

gene=""
print(g.readline()) 
for line in g:
    gene=gene+line[:-1]
    
gene


# In[58]:


r = open("rna.fna", "r")

rna=""
print(r.readline()) 
for line in r:
    rna=rna+line[:-1]
    
rna


# In[16]:





# let's find length of each exon and therefore each intron
# 
# using biopython library to visualize gene https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec4
# 
# chapter 17 talks about graphic

# In[22]:





# In[ ]:




