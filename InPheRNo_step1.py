# -*- coding: utf-8 -*-
"""
@author: emad2

This is the first step of InPheRNo. This script uses the normalized FPKM data for TCGA and generats P-value 
output files necessary to use in the next steps. In calculating phenotype-gene
p-values, we are not using any gene-gene network information. The gene-TF p-values
are calculated either using Pearson correlation, or using ElasticNet.

As input, this script take in FPKM values in the form of a gene x samples matrix. 
The outputs are placed in a folder called PGM_data. If tf_gene_method == ElasticNet,
the matrix of TF-gene p-values assign a value -1 to entries that were not selected
based on Elastic Net. Note that the gene-pheno p-values are sorted such that smallest 
p-vals appear first. Also TF-gene matrix is sorted so that gene names match the gene-pheno
matrix.

"""

import scipy.special as ssp
import numpy as np
import glob
import pandas as pd
import time
import os
import scipy.stats as ss
from sklearn.linear_model import ElasticNetCV
import statsmodels.api as sm
import pickle

address_TF = '/Users/emad2/Amin_Research/TRN/Data/TF_Ensemb_AnimalTFDB.txt'

address_out_dir = '/Users/emad2/Amin_Research/TRN/Data/TCGA_DataCommons/FPKM_Data/PGM_data'
address_in_expr_GTEX = '/Users/emad2/Amin_Research/TRN/Data/GTEX_V6p/GTEX_EXPR_RPKM_IQN.csv'
address_in_expr_TCGA = '/Users/emad2/Amin_Research/TRN/Data/TCGA_DataCommons/FPKM_Data/TCGA_EXPR_FPKM_IQN.csv'



tf_gene_method = 'ElasticNet'  #Options are Pearson and ElasticNet
if tf_gene_method == 'ElasticNet':
    max_num_coefs = 15  #since we want to limit number of nonzero coefs, this max value is useful
    l1_rat = 0.5    #l1_ratio used for ElasticNet

print(tf_gene_method)

eps = 3e-308


TF_list = list(pd.read_csv(address_TF, sep='\t', header=None, index_col=0).index.values)

address_out_pvalue_gt_generic = os.path.join(address_out_dir, 'Pvalue_gene_TF_%s_generic_FPKM.csv' %(tf_gene_method))

address_pickle_out = os.path.join(address_out_dir, 'Pvalue_gene_TF_generic_%s.pickle' %(tf_gene_method))
#I may not want to include AdrenalGland_ACC
cancer_list = ['AdrenalGland_ACC', 'AdrenalGland_PCPG', 'Brain_GBM', 
               'Brain_LGG', 'Breast', 'Colorectal_COAD', 'Colorectal_READ',
               'Esophagus', 'Liver', 'Lung_LUAD', 'Lung_LUSC', 'Ovary',
               'Pancreas', 'Prostate', 'Skin', 'Stomach', 'Testis', 'Thyroid']
#cancer_list = ['Breast', 'Ovary']
#cancer_list = ['AdrenalGland_ACC', 'AdrenalGland_PCPG']





###############################################################################
sample_name = []    #a list of lists where each list contains samples in a cancer in cancer_list
expr_in = pd.DataFrame()
list_ = []


for ic in range(len(cancer_list)):
    cancer = cancer_list[ic]
    print(cancer)
    address_in_dir = '/Users/emad2/Amin_Research/TRN/Data/TCGA_DataCommons/FPKM_Data/%s_GDC/Processed' %cancer
    address_out_dir = '/Users/emad2/Amin_Research/TRN/Data/TCGA_DataCommons/FPKM_Data/PGM_data'
    address_in_FPKM_primary = os.path.join(address_in_dir, 'FPKM_PrimaryTumor_%s.csv' %cancer)
    
    df = pd.read_csv(address_in_FPKM_primary, index_col=0, header=0)    #reads the FPKM values for each cancer
    list_.append(df)    #list_ contains dataframes of different cancers
    sample_name.append(list(df.columns))


##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
expr_IQN_GTEX = pd.read_csv(address_in_expr_GTEX, index_col=0)    #this is IQN of log2(RPKM+1) values for all cancers (not sorted)
expr_IQN = pd.read_csv(address_in_expr_TCGA, index_col=0)    #this is IQN of log2(RPKM+1) values for all cancers (not sorted)
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########

gene_TF_list = list(expr_IQN.index.values)
gene_TF_list_GTEX = list(expr_IQN_GTEX.index.values)
gene_TF_list = list(set(gene_TF_list).intersection(gene_TF_list_GTEX))


only_TF_list = list(set(TF_list).intersection(set(gene_TF_list)))
only_TF_list.sort()
only_gene_list = list(set(gene_TF_list) - set(TF_list))
only_gene_list.sort()

expr_IQN_gene = expr_IQN.loc[only_gene_list]    #this contains expression of all samples
expr_IQN_TF = expr_IQN.loc[only_TF_list]


#Form a dataframe of gene x TF for pvalue_gt. This DF will be row-sorted depending on the cancer
pvalue_gt_array = (-1) * np.ones((len(only_gene_list), len(only_TF_list))) #A gene x TF matrix


if tf_gene_method == 'Pearson':
    for i in range(len(only_gene_list)):
        print('Pvalue_gene_TF', i)
        for j in range(len(only_TF_list)):
            _, pvalue_gt_array[i][j] = ss.pearsonr(expr_IQN_gene.iloc[i], expr_IQN_TF.iloc[j])
            if pvalue_gt_array[i][j] < eps:
                pvalue_gt_array[i][j] = eps
elif tf_gene_method == 'ElasticNet':
    X_features = expr_IQN_TF.values.T
    start_time = time.clock()

    for i in range(len(only_gene_list)):
        print('Pvalue_gene_TF', i)
        y = expr_IQN_gene.iloc[i].values
        EN_model = ElasticNetCV(l1_ratio=l1_rat)

        ####make sure that number of nonzero coefs do not exceed max_num_coefs
        alphas1, coefs1, _ = EN_model.path(X_features, y, eps=0.01, n_alphas=10)
        num_coefs = np.sum(coefs1 != 0, axis=0)
        #print(num_coefs)
        #print(num_coefs[num_coefs <= max_num_coefs][-1])
        rep_EN = 0
        while (num_coefs[0] != num_coefs[-1]) and (num_coefs[num_coefs <= max_num_coefs][-1] != max_num_coefs) and (rep_EN < 10):
            rep_EN += 1
            alpha_min = alphas1[(num_coefs <= max_num_coefs)][-1]
            alpha_max = alphas1[(num_coefs > max_num_coefs)][0]
            alphas3 = np.linspace(alpha_min, alpha_max, 10)
            alphas1, coefs1, _ = EN_model.path(X_features, y, alphas=alphas3)
            num_coefs = np.sum(coefs1 != 0, axis=0)
            #print(num_coefs)
            #print(num_coefs[num_coefs <= max_num_coefs][-1])        
        print('repeat', rep_EN )
        if num_coefs[0] == num_coefs[-1]:
            EN_coef = coefs1[:, 0]
            selected_ind = np.array(range(len(only_TF_list)))[EN_coef!=0]   
        else:
            EN_coef = coefs1[:, len(num_coefs[num_coefs <= max_num_coefs])-1]
            selected_ind = np.array(range(len(only_TF_list)))[EN_coef!=0]
            
        print(len(selected_ind))
        n_selected_TF = len(selected_ind)
        X_features_new = X_features[:, selected_ind]
     
        model = sm.OLS(y, X_features_new)
        results = model.fit()
     
        ts_b = results.tvalues
        

        #####prcise estimation of p-vals
        N = np.shape(X_features)[0]-1   #degrees of freedom in ttest for regression
        pvals_precise = []
        x1_bar = 0  #mean of the first distribution
        n1 = (N+2)//2   #num obs first distribution
        n2 = N+2-n1     #num obs second distribution
        x2_std = 0
        t_precise = []
        for j in range(len(ts_b)):
            T = ts_b[j]     #T statistic
            x2_bar = -np.sign(T)    #mean of the second distribution
            x1_std=1/T*np.sqrt((n1+n2-2)/(n1-1)/(1/n1+1/n2))    #std of first distribution
            (t, p) = ss.ttest_ind_from_stats(x1_bar, x1_std, n1, x2_bar, x2_std, n2)
            if p < eps:
                p = eps
            pvals_precise.append(p)
            t_precise.append(t)
        
        pvalue_gt_array[i, selected_ind] = pvals_precise
        with open(address_pickle_out, "wb") as f:
            pickle.dump([pvalue_gt_array, i], f)

pvalue_gt_df = pd.DataFrame(pvalue_gt_array, index=only_gene_list, columns=only_TF_list)

##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
pvalue_gt_df.to_csv(address_out_pvalue_gt_generic) #Not sorted gene-TF correlation p-values
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########


for ci in range(len(cancer_list)):
    cancer = cancer_list[ci]
    pvalue_gp_array = np.zeros(len(only_gene_list))
    sample_interest = sample_name[ci]
    sample_other = []
    for cj in range(len(cancer_list)):
        if cj != ci:
            sample_other += sample_name[cj] #list of all samples that are not in cancer of interest
    for i in range(len(only_gene_list)):
        _, pvalue_gp_array[i] = ss.ttest_ind(expr_IQN_gene.loc[only_gene_list[i]][sample_interest], expr_IQN_gene.loc[only_gene_list[i]][sample_other], equal_var=False)
        if pvalue_gp_array[i] < eps:
            pvalue_gp_array[i] = eps
    print(cancer, '# DEG %s' %(pvalue_gp_array < 0.05).sum())

    #sort all genes based on pvalues and form dataframe:
    argsort_ind = np.argsort(pvalue_gp_array)
    pvalue_gp_array_psorted = pvalue_gp_array[argsort_ind]
    only_gene_list_psorted = [only_gene_list[i] for i in argsort_ind]
    pvalue_gp_df_sorted = pd.DataFrame(pvalue_gp_array_psorted, index=only_gene_list_psorted, columns=['PValue'])   #this is the dataframe to be written to file

    #Now we form gene_TF dataframe where genes are sorted based on only_gene_list_psorted
    pvalue_gt_df_sorted = pvalue_gt_df.loc[only_gene_list_psorted]
    
    #write to file:
    address_pval_gp = os.path.join(address_out_dir, 'Pvalue_gene_phenotype_1vsAll_FPKM_%s.csv' %(cancer))
    address_pval_tg = os.path.join(address_out_dir, 'Pvalue_TF_gene_%s_1vsAll_FPKM_%s.csv' %(tf_gene_method, cancer))
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
    pvalue_gp_df_sorted.to_csv(address_pval_gp) #sorted pvalues of gene-pheno
    pvalue_gt_df_sorted.to_csv(address_pval_tg) #sorted pvalues of gene-TF (sorted based on gene-pheno)
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########


