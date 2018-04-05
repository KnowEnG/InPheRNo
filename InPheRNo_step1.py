# -*- coding: utf-8 -*-
"""
@author: emad2

This is the first step of InPheRNo. This script uses the normalized transcriptomic
data to generate P-value of gene-TF associations using ElasticNet. The results
generated for each (TF,gene) pair are pseudo-pvalues that approximately follow
a Beta distribution under the Null hypothesis. 


As input, this script takes in 3 files: 1) A list of ALL TFs, 2) a file containing
p-values of gene-phenotype associations only for genes of interest, and 3) the 
expression profiles of genes and TFs. If there are shared gene name sbetween list
of TFs and gene-phenotype pvalue file, the script drops those. 


The script generates two outputs: a gene-phenotype p-value file and a gene-TF
pseudo p-value file such that the gene names match in the two files. These files
only include genes of interest and not all genes. The matrix of TF-gene p-values
contains values of -1 for entries which were assigned a coefficient of 0 in 
ElasticNet.


"""

#import scipy.special as ssp
import numpy as np
#import glob
import pandas as pd
import time
import os
import scipy.stats as ss
from sklearn.linear_model import ElasticNetCV
import statsmodels.api as sm
#import pickle
import argparse

###############################################################################
# Parse command line options/arguments
parser = argparse.ArgumentParser()

parser.add_argument('-id', '--input_directory', default = './Data', help = 'Address of directory containing input files')
parser.add_argument('-od', '--output_directory', default = './Results', help = 'output directory adddress')
parser.add_argument('-it', '--input_tf', default = 'TF_Ensemble.csv', help = 'Name of the file containing list of TFs in a csv file. The file should not have a header.')
parser.add_argument('-ie', '--input_expression', default = 'expr_sample.csv', help = 'A file containing gene and TF expression data (gene x samples). The file has a header (sample names).')
parser.add_argument('-igp', '--input_gene_phenotype_interest', default = 'Pvalue_gene_phenotype_interest.csv', help = 'A file (gene x pvalue) containing p-values of gene-phenotype only for genes of interest (and not all genes), sorted in an ascending order based on the p-value (smallest p-values appear first). Only include genes of interest to reduce computation time. The file has a header.')
parser.add_argument('-mt', '--max_num_tf', default = 15, help = 'Maximum number of Tfs recovered for each gene using EN')
parser.add_argument('-lr', '--l1_ratio', default = 0.5, help = 'l1 ratio in EN model' )
parser.add_argument('-ogp', '--output_gene_phenotype', default = 'Pvalue_gene_phenotype_interest_tmp.csv', help = 'A file (gene x pvalue) containing p-values of gene-phenotype, sorted in an ascending order based on the p-value (smallest p-values appear first). The file has a header.')
parser.add_argument('-tgt', '--output_gene_tf', default = 'Pvalue_gene_tf_tmp.csv', help = 'A file (gene x tf) containing p-values of gene-tf, sorted in an ascending order based on gene-pheno file. The file has a header.')

args = parser.parse_args()

###############################################################################
delim_tl = ','
delim_ex = ','
delim_gp = ','
if args.input_tf[-3:] in ['tsv', 'txt']:
    delim_tl = '\t'
if args.input_expression[-3:] in ['tsv', 'txt']:
    delim_ex = '\t'
if args.input_gene_phenotype_interest[-3:] in ['tsv', 'txt']:
    delim_gp = '\t'



address_TF = os.path.join(args.input_directory, args.input_tf)
address_out_dir = args.output_directory
if not os.path.exists(address_out_dir):
    os.makedirs(address_out_dir)

address_in_expr = os.path.join(args.input_directory, args.input_expression) 
address_in_gene_pheno = os.path.join(args.input_directory, args.input_gene_phenotype_interest) 

address_out_pvalue_gp= os.path.join(args.output_directory, args.output_gene_phenotype)
address_out_pvalue_gt = os.path.join(args.output_directory, args.output_gene_tf)



max_num_coefs = args.max_num_tf  #since we want to limit number of nonzero coefs, this max value is useful
l1_rat = args.l1_ratio    #l1_ratio used for ElasticNet

np.random.seed(1011)
eps = 3e-308



##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
TF_list = list(pd.read_csv(address_TF, sep=delim_tl, header=None, index_col=0).index.values)
expr_all = pd.read_csv(address_in_expr, sep=delim_ex, index_col=0)    #this is IQN of log2(RPKM+1) values for all cancers (not sorted)
gene_pheno_pval = pd.read_csv(address_in_gene_pheno, sep=delim_gp, index_col=0)
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########


gene_pheno_list = list(gene_pheno_pval.index.values)
gene_TF_list = list(expr_all.index.values)
#gene_TF_list_GTEX = list(expr_IQN_GTEX.index.values)
#gene_TF_list = list(set(gene_TF_list).intersection(gene_TF_list_GTEX))


only_TF_list = list(set(TF_list).intersection(set(gene_TF_list)))
only_TF_list.sort()
only_gene_list = list(set(gene_TF_list) - set(TF_list))
only_gene_list = list(set(only_gene_list).intersection(gene_pheno_list))
only_gene_list.sort()


#find intersetion of gene_pheno_pval file and genes with expression values
gene_pheno_pval = gene_pheno_pval.loc[only_gene_list]
pvalue_gp_df = gene_pheno_pval.sort_values(by=gene_pheno_pval.columns[0])    #sort
only_gene_list = list(pvalue_gp_df.index.values)



expr_gene = expr_all.loc[only_gene_list]    #this contains expression of all samples
expr_TF = expr_all.loc[only_TF_list]


#Form a dataframe of gene x TF for pvalue_gt. This DF will be row-sorted depending on the cancer
pvalue_gt_array = (-1) * np.ones((len(only_gene_list), len(only_TF_list))) #A gene x TF matrix


X_features = expr_TF.values.T
start_time = time.clock()

for i in range(len(only_gene_list)):
    print('Pvalue_gene_TF', i)
    y = expr_gene.iloc[i].values
    EN_model = ElasticNetCV(l1_ratio=l1_rat)

    ####make sure that number of nonzero coefs do not exceed max_num_coefs
    alphas1, coefs1, _ = EN_model.path(X_features, y, eps=0.01, n_alphas=10)
    num_coefs = np.sum(coefs1 != 0, axis=0)
    #print(num_coefs)
    #print(num_coefs[num_coefs <= max_num_coefs][-1])
    rep_EN = 0
    #    while (num_coefs[0] != num_coefs[-1]) and (num_coefs[num_coefs <= max_num_coefs][-1] != max_num_coefs) and (rep_EN < 10):
    while (num_coefs[0] != num_coefs[-1]) and (max(num_coefs[num_coefs <= max_num_coefs]) != max_num_coefs) and (rep_EN < 10):
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


pvalue_gt_df = pd.DataFrame(pvalue_gt_array, index=only_gene_list, columns=only_TF_list)




##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
pvalue_gt_df.to_csv(address_out_pvalue_gt) #sorted gene-TF pseudo p-values
pvalue_gp_df.to_csv(address_out_pvalue_gp)
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########


