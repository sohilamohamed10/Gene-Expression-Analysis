# Gene-Expression-Analysis
## Analysis of gene expression (GE) data for the cancer type Lung Squamous Cell Carcinoma (LUSC)

## Steps 
### a)Correlation
Computing the correlation between the normal samples and the diseased samples for
each gene.
### Results 
Correlation coeffients.

 ![](Results/correlationCoeff.png)

Ranked genes based on their correlation coefficient and the highest positive CC and the lowest negative CC and the names of these two genes.

 ![](Results/RankedCC.png) 

 Plotting gene expressions for these 2 genes.
  
 ![](Results/correlationPlot1.png) 
 
 ![](Results/correlationPlot2.png) 


### a)Hypothesis testing
Infer the differentially expressed genes (DEGs); the genes whose
expression level differ from one condition (healthy) to another (diseased).
### Results
Applying t_test and computing  p_values and level of significance for two cases paired and independent samples. 

Applying FDR correction for both cases.

 ![](Results/pairedP_values.png) 

 ![](Results/indP_values.png) 


 Compare the two DEGs sets (paired and independent) after the FDR correction common ,distinct paired and distinct independent 

 ![](Results/common&distinctPairedP_values.png) 

 ![](Results/distinctIndP_values.png) 


