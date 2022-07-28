#!/usr/bin/env python
# coding: utf-8

# In[42]:


import pandas as pd
import json

#############################################
''' 
    get_orientation(geneName)
        Extracts the orientation of the sequence from 
        the .tsv report for the given gene
        params: gene name(str)
        output: one of the following strings: +;+ , -;- , -
'''
def get_orientation(geneName):
    table = pd.read_csv(geneName+ "/data_table.tsv", sep='\t')
    return table["orientation"][0]


#############################################
'''
    annotate_gene(geneName,accessVer)
        Given the name of the gene and desired accession version,
        produce a dict w/ position info for each exon and intron
        params: gene name(str), accession version(str)
        output: a dictionary 
                key=exon/intron + #
                value=list(length, pos, pos on chr, pos on rna-seq)
'''
def annotate_gene(geneName,accessVer):
    
    re = open(geneName+ "/data_report.jsonl")
    report = json.load(re)
    
    for i in range(len(report['transcripts'])): 
        if report['transcripts'][i]['accessionVersion'] ==  accessVer:
            exons = report['transcripts'][i]['genomicLocations'][0]['exons'] ##extracting exon position of the desired accession version 
           
    ##fetch the orientation of sequence
    orientation = get_orientation(geneName)
    ann_gene = dict() #the outout
                      #key=exon/intron # : value=list(length, pos, pos on chr, pos on rna-seq)
    count=1  #counts the number of exons
    s=1  #starting position
    rna_begin=1
    
    for i in range(len(exons)):
        
        if orientation == "+;+":
            
            #exon position
            begin=int(exons[i]["begin"])
            end=int(exons[i]["end"])
            lenEx= end - begin  + 1

            #adding this exon plus its begin and end pos to the dict
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

        elif orientation in ["-", "-;-"]:
            
            #exon position
            begin=int(exons[i]["end"])
            end=int(exons[i]["begin"])
            lenEx = begin - end  + 1
            
            #adding this exon plus its begin and end pos to the dict
            ann_gene["exon " + str(count)]=[lenEx,(s,s+lenEx-1),(begin,end),(rna_begin, rna_begin+lenEx-1)]
            rna_begin=rna_begin+lenEx

            #intron positions
            if (count < len(exons)):
                #intron
                intEnd=int(exons[i+1]["end"])+1
                intBegin=end-1

                lenInt= intBegin - intEnd + 1

                ann_gene["intron " + str(count)]=[lenInt,(s+lenEx,s+lenEx+lenInt-1),(intBegin,intEnd),"-"]

                s=s+lenEx+lenInt
                count=count+1
    
    return ann_gene

#############################################
'''
    make_df_table(gene_dict)
        Given the dictionary output of annotate_gene(), or a dict of
        same format, it will transform it into pandas dataframe
        params:gene_dict(dict)
        output:pandas df
'''
def make_df_table(gene_dict):
    gene_df=pd.DataFrame(gene_dict).T
    gene_df.columns=["exon/intron len", "position", "pos on chrom", 
                "pos on rna-seq"]
    return gene_df



if __name__=="__main__":
    gene=annotate_gene("tyrobp","NM_003332.4")
    tb=make_df_table(gene)
    print(tb)

