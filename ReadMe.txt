This folder contains the R code and data for paper "Truncated Rank-Based Tests for Two-Part Models with Excessive Zeros and Applications to Microbiome Data" by Wang, Chen and Li. 

It contains the 8 code files and two data files. 

#######################################################
Code files: 

ARE_Wilcoxon.R: it is to generate the theoretical ARE between the truncated Wilcoxon statistics and standard Wilcoxon statistics. This file will generate Figure 2 in the paper. 

ARE_ANOVA.R: it is to generate the theoretical ARE between the truncated Kruskal-Wallis statistics and standard Kruskal-Wallis statistics. This file will generate Figure 3 in the paper. 

Simulation.R: This file executes the simulation experiments. It will generate 7 Rdata files that contains the simulated p-values and parameters. It will also generate Figure 4 in the paper. 

Normal_vs_disease.R: This file is to do the comparison of gut microbiome between normal and Crohn’s disease patients. 
It requires "Wilcoxon.zeros.R", and the data files "G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Disease" and "G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Normal".
It will generate Figure 5 in the paper. 

AntiTNF_Different_Time.R: This file is to do the comparison of gut microbiome across time after treatment. 
It requires "ANOVA.zeros.R", and the data files "G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Disease"  and "COMBO_PLEASE_Sample_Information_Taxa".
It will generate Figure 6 in the paper. 

AREsimulation.R: This file will simulate the ARE curve, and compare with the theoretical results. It will create Figure 1 in our supplementary materials.


Wilcoxon.zeros.R: the function of our truncated statistic.

ANOVA.zeros.R: the function of our truncated statistic. 

#######################################################
Data files: 

G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Disease.xls

G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Normal.xls

COMBO_PLEASE_Sample_Information_Taxa.csv

