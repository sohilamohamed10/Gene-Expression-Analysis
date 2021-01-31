#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd 
from pandas import DataFrame
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind,ttest_rel
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests


# In[ ]:


#loading files & filtration
import pandas as pd
healthy=pd.read_csv('healthy.txt', sep='\t')
cancer=pd.read_csv('cancer.txt', sep='\t')
filtered_h=healthy[healthy.astype('bool').mean(axis=1)>=0.5] 
filtered_c=cancer[cancer.astype('bool').mean(axis=1)>=0.5]
common = filtered_h.merge(filtered_c,on=['Hugo_Symbol'])
final_healthy = common.iloc[:, :52]
final_cancer = common.iloc[:, 52:]
final_cancer.insert(0, 'Hugo_Symbol',final_healthy.iloc[0:,0]) 


# In[ ]:


# computing correlation coefficients and listing them descendingly
mylist=[]
gene_idx=[]
for i in range (0 , len(final_healthy)):
	gene_h=final_healthy.iloc[i,2:]
	gene_c=final_cancer.iloc[i,2:]
	r, _= pearsonr(gene_h,gene_c)
	mylist.append(r)
	gene_idx.append(i)
   
list_of_tuples = list(zip(mylist,gene_idx))  
df=pd.DataFrame(list_of_tuples,columns=["correlation coeffs.","gene index"])
sorted_cc=df.sort_values('correlation coeffs.',ascending=False)
print(sorted_cc)
mylist.sort(reverse=True)
ranked_cc=pd.DataFrame(mylist,columns=["correlation coeffs."])
print(ranked_cc) 

# In[ ]:


#printing genes of highest and lowest correlation coefficient
h_cc_idx=sorted_cc.iloc[0,1]
print(h_cc_idx)
l_cc_idx=sorted_cc.iloc[len(sorted_cc)-1,1]
highest_cc=sorted_cc.iloc[0,0]
highest_gene=final_healthy.iloc[h_cc_idx,0]
lowest_cc=sorted_cc.iloc[len(sorted_cc)-1,0]
lowest_gene=final_healthy.iloc[l_cc_idx,0]
correlation=pd.DataFrame({'gene':[highest_gene,lowest_gene],'CC':[highest_cc,lowest_cc]})
correlation


# In[ ]:


#ploting gene expression of highest cc gene
hh_cc=final_healthy.iloc[h_cc_idx,2:]
hc_cc=final_cancer.iloc[h_cc_idx,2:]
data_h = list(zip(hh_cc,hc_cc)) 
print(hh_cc) 
data_high=pd.DataFrame(data_h,columns=["hh_cc","hc_cc"])
data_high.plot(x ='hh_cc', y='hc_cc', kind = 'scatter')
plt.show()

# In[ ]:


#ploting gene expression of lowest cc gene
lh_cc=final_healthy.iloc[l_cc_idx,2:]
lc_cc=final_cancer.iloc[l_cc_idx,2:]
data_l = list(zip(lh_cc,lc_cc))  
data_low=pd.DataFrame(data_l,columns=["lh_cc","lc_cc"])
data_low.plot(x ='lh_cc', y='lc_cc', kind = 'scatter')
plt.show()



# In[ ]:


# Hypothesis testing
Genes_name=final_healthy["Hugo_Symbol"]  #get the gene_name from the file
Genes=[] # list of the genes_name
Related_PValue=[] #list of the calculated paired_pvalue for each gene
Independent_PValue=[]  #list of the calculated independent_pvalue for each gene
significance_RPVal=[]  #list of the significance of the paired pvalue 
significance_INDPVal=[] #list of the significance of the independent pvalue
for i in range(len(final_healthy)): #for loop to iterate on the rows
    Gi_healthy=final_healthy.iloc[i,2:]  #get the G_i for all columns from the healtyh samples
    Gi_cancer=final_cancer.iloc[i,2:]    #get the G_i for all columns from the cancer samples
    RP_Value=ttest_rel(Gi_healthy,Gi_cancer).pvalue #calculate the paired_pvalue for each gene
    IndP_Value=ttest_ind(Gi_healthy,Gi_cancer).pvalue #calculate the independent_pvalue for each gene
    Genes.append('G'+str(i))  
    Related_PValue.append(RP_Value)
    Independent_PValue.append(IndP_Value)
    if RP_Value < 0.05:
        significance_RPVal.append("True")
    else:
        significance_RPVal.append("False")
    if IndP_Value < 0.05:
        significance_INDPVal.append("True")
    else:
        significance_INDPVal.append("False")

Paired_PV=pd.DataFrame({'Gene_name':Genes_name,'Rel_PValue':Related_PValue,'significance_RPVal':significance_RPVal})
Independent_PV=pd.DataFrame({'Gene_name':Genes_name,'Ind_PValue':Independent_PValue,'significance_INDPVal':significance_INDPVal})
DEG_p_values_rel_before = Paired_PV.loc[Paired_PV['significance_RPVal']== 'True']
DEG_p_values_rel_before
#DEG_p_values_ind_before = Independent_PV.loc[Independent_PV['significance_INDPVal']== 'True']
#DEG_p_values_ind_before
#Paired_PV
#Independent_PV 


# In[ ]:


#correction THE PAIRED by FDR
corrected_p_values_rel = multipletests(Related_PValue, alpha=0.05, method='fdr_bh')[1]
correct_significance_rel = multipletests(Related_PValue, alpha=0.05, method='fdr_bh')[0]
all_corrected_p_values_rel=pd.DataFrame({'Gene_name':final_healthy.iloc[0:,0], 'corrected_p_values_rel ':corrected_p_values_rel ,'correct_significance_rel':correct_significance_rel}) 
DEG_corrected_p_values_rel= all_corrected_p_values_rel[all_corrected_p_values_rel['correct_significance_rel']== True]

print(all_corrected_p_values_rel)


# In[ ]:


#correction THE INDEPENDENT by FDR
corrected_p_values_ind = multipletests(Independent_PValue, alpha=0.05, method='fdr_bh')[1]
correct_significance_ind = multipletests(Independent_PValue, alpha=0.05, method='fdr_bh')[0]
all_corrected_p_values_ind=pd.DataFrame({'Gene_name':final_healthy.iloc[0:,0], 'corrected_p_values_ind ':corrected_p_values_ind ,'correct_significance_ind':correct_significance_ind}) 
DEG_corrected_p_values_ind= all_corrected_p_values_ind[all_corrected_p_values_ind['correct_significance_ind']== True]

print(all_corrected_p_values_ind)


# In[ ]:


#COMMON GENES between fdr_indpendent and fdr_paired
common_DEG = DEG_corrected_p_values_ind.merge(DEG_corrected_p_values_rel,on=['Gene_name'])
#common_DEG
gene_common_DEG=common_DEG.iloc[0:,0]
gene_common_DEG


# In[ ]:


#distinct fdr_paird
distinct_rel = DEG_corrected_p_values_rel.merge(common_DEG, how = 'outer' ,indicator=True).loc[lambda x : x['_merge']=='left_only']
#distinct_rel['Gene_name'].to_list()
DEG_distinct_rel=distinct_rel.iloc[0:,0]
DEG_distinct_rel


# In[ ]:


#distinct fdr_INDPENDENT
distinct_ind = DEG_corrected_p_values_ind.merge(common_DEG, how = 'outer' ,indicator=True).loc[lambda x : x['_merge']=='left_only']
#distinct_ind
#distinct_ind['Gene_name'].to_list()
gene_distinct_ind=distinct_ind.iloc[0:,0]
gene_distinct_ind


# In[ ]:




