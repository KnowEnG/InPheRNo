# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:03:49 2017

@author: emad2

This is the first step in a series of steps in the phenotype-relevant TRN 
recovery project. This script uses the normalized FPKM data for TCGA and generats P-value 
output iles necessary to use in the next steps. In calculating phenotype-gene
p-values, we are using any gene-gene network information. The gene-TF p-values
are calculated either using Pearson correlation, or using ElasticNet.

As input, this script take in FPKM values in the form of a gene x samples matrix. 
The outputs are placed in a folder called PGM_data. If tf_gene_method == ElasticNet,
the matrix of TF-gene p-values assign a value -1 to entries that were not selected
based on Elastic Net. Note that the gene-pheno p-values are sorted such that smallest 
p-vals appear first. Also TF-gene matrix is sorted so that gene names match the gene-pheno
matrix.

"""

import numpy as np
import pandas as pd
import os
import scipy.stats as ss
import pickle
from scipy import sparse
from sklearn.preprocessing import normalize



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

#I may not want to include AdrenalGland_ACC
cancer_list = ['AdrenalGland_ACC', 'AdrenalGland_PCPG', 'Brain_GBM', 
               'Brain_LGG', 'Breast', 'Colorectal_COAD', 'Colorectal_READ',
               'Esophagus', 'Liver', 'Lung_LUAD', 'Lung_LUSC', 'Ovary',
               'Pancreas', 'Prostate', 'Skin', 'Stomach', 'Testis', 'Thyroid']
#cancer_list = ['Breast', 'Ovary']
#cancer_list = ['AdrenalGland_ACC', 'AdrenalGland_PCPG']




##############################################################################
###Network info
network_name = 'STRING_experimental'
restart_prob = 0.5  # Restart probability in RWR
tolerance = 1e-8    # Tolerance for RWR
max_iter = 100      # Maximum number of iterations in RWR

address_net = os.path.join('/Users/emad2/Amin_Research/Drug_Response/KN', network_name + '.csv')



###############################################################################
def is_number(num):
    """
    Determine whether a string s is a number (i.e., any floating point
    representation of a number, including scientific notation)
    """
    try:
        float(num)
        return True
    except ValueError:
        return False



###############################################################################
def SmatchN(expr_DF_in, node_names_in):
    """
    Matches S (spreadsheet of gene expressions) and N (network)
    The function returns expr_DF_out which is formed by reshuffling columns of
    expr_DF_in. Also, node_names_out is formed by reshuffling node_names_in. The
    intersection of node_names_out and column names of expr_DF_out are placed at
    the beginning of both lists.
    """
    node_names_in_set = set(node_names_in)
    gene_names_in_set = set(expr_DF_in.columns.values)

    nodes_genes_intersect = sorted(list(gene_names_in_set & node_names_in_set))
    nodes_minus_genes = sorted(list(node_names_in_set - gene_names_in_set))
    genes_minus_nodes = sorted(list(gene_names_in_set - node_names_in_set))

    genes_names_out = nodes_genes_intersect + genes_minus_nodes
    nodes_names_out = nodes_genes_intersect + nodes_minus_genes
    expr_DF_out = expr_DF_in [genes_names_out]
    return(expr_DF_out, nodes_names_out, nodes_genes_intersect)


###############################################################################
def RWR_matrix(node_names, network_matrix, restart_matrix, restart_prob, max_iter, tolerance):
    """Performs a RWR (Random Walk with Restart) with the given parameters"""
    no_restart_prob = 1 - restart_prob
    init_prob = 1/len(node_names)

    # Create the vector of probabilities for the nodes
    steady_prob_old = np.empty(np.shape(restart_matrix))
    steady_prob_old.fill(init_prob)

    residual = 100
    num_iter = 0
    while (residual > tolerance) and (num_iter < max_iter):
        steady_prob_new = (sparse.csr_matrix.dot(steady_prob_old, network_matrix) * no_restart_prob
                            + restart_prob * restart_matrix)

        residual = max(abs(steady_prob_new - steady_prob_old).sum(axis=1))
        num_iter += 1
        print('RWR_iteration = ', num_iter)

        steady_prob_old = steady_prob_new.copy()
    return(num_iter, residual, steady_prob_new)




##############################################################################
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
###This is where commentable section starts
"""
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



###############################################################################
DEFAULT_COLUMN_HEADERS = ['n_alias_1', 'n_alias_2', 'weight', 'type']
delimiter = ','
# Step 1: Read the input
# the input_file is the network
input_file = address_net
with open(input_file, 'r') as f:
    # Check whether the first line is data or headers
    first_line = f.readline().strip()
    f.seek(0)   #go back to the beginning of the file
    fields = first_line.split(sep=delimiter)
    # data
    if is_number(fields[2]):
        net_DF = pd.read_csv(f, sep=delimiter, names=DEFAULT_COLUMN_HEADERS)
    # headers
    else:
        net_DF = pd.read_csv(f, sep=delimiter, header=0)

#Remove TF
net_DF = net_DF.loc[~net_DF['n_alias'].isin(set(TF_list)) & ~ net_DF['n_alias.1'].isin(set(TF_list))]   #keep only genes and exclude TFs from network


# Get the column headers
n1 = net_DF.columns[0]
n2 = net_DF.columns[1]
w = net_DF.columns[2]
#        t = net_DF.columns[3]

# Get the unique nodes -- the first two columns of the input data,
# converted to sets to remove duplicates, union'ed, then sorted
nodes1 = net_DF.iloc[:, 0]
nodes2 = net_DF.iloc[:, 1]
nodes = set(nodes1) | set(nodes2)
node_names = sorted(nodes)
num_nodes = len(node_names)
#print(node_names)

# Output some info about the input data
print("Input:")
print("Lines of data:", len(net_DF))
print("Number of unique nodes:", num_nodes)
print()




##############################################################################

expr_IQN = pd.read_csv(address_in_expr_TCGA, index_col=0)    #this is IQN of log2(RPKM+1) values for all cancers (not sorted)



gene_TF_list = list(expr_IQN.index.values)

only_TF_list = list(set(TF_list).intersection(set(gene_TF_list)))
only_TF_list.sort()
only_gene_list = list(set(gene_TF_list) - set(TF_list))
only_gene_list.sort()

expr_IQN_gene = expr_IQN.loc[only_gene_list]    #this contains expression of all samples
expr_IQN_TF = expr_IQN.loc[only_TF_list]


##############################################################################
###Mapping
# Reorder gene column names of expression spreadsheet and gene node names of network
# such that both name lists start with the intersection of the two list ordered alphabetically
(expr_IQN_gene, node_names, nodes_genes_intersect) = SmatchN(expr_IQN_gene.T, node_names)
expr_IQN_gene = expr_IQN_gene.T
node2index = {node_names[i]:i for i in range(len(node_names))}
index2node = node_names.copy()

only_gene_list = list(expr_IQN_gene.index.values)

gene2index = {only_gene_list[i]:i for i in range(len(only_gene_list))}
index2gene = only_gene_list.copy()





###############################################################################
# Transform the first two columns of the DataFrame -- the nodes -- to their indexes
net_DF[n1] = net_DF[n1].apply(lambda x: node2index[x])
net_DF[n2] = net_DF[n2].apply(lambda x: node2index[x])

# Create the sparse matrix
network_matrix = sparse.csr_matrix((net_DF[w].values, (net_DF[n1].values, net_DF[n2].values)),
                                   shape=(num_nodes, num_nodes), dtype=float)
# Make the ajdacency matrix symmetric
network_matrix = (network_matrix + network_matrix.T)
network_matrix.setdiag(0)

# Normalize the rows of network_matrix because we are multiplying vector by matrix (from left)
network_matrix = normalize(network_matrix, norm='l1', axis=1)

restart_matrix = np.eye(len(node_names),len(nodes_genes_intersect)).T


# For the genes if interest (intersection of genes in network and spreadsheet),
# this RWR produces similarity with other genes
(num_iter, residual, gene_similarity_smooth) = RWR_matrix(node_names, network_matrix, restart_matrix, restart_prob, max_iter, tolerance)

gene_similarity_smooth = gene_similarity_smooth[:, 0:np.size(gene_similarity_smooth, axis=0)]
gene_similarity_smooth = normalize(gene_similarity_smooth, norm='l1', axis=1)





###############################################################################
###############################################################################
### TEMP: this is written simply to save the network smoothed expression as a pickle so that 
### we do not have to repeat the process many times.

net_smooth_pickle = os.path.join(address_out_dir, 'Similarity_TCGA_RWR%s_%s' %(int(100 * restart_prob), network_name))

with open(net_smooth_pickle,"wb") as f:
    pickle.dump([nodes_genes_intersect, sample_name, expr_IQN_gene, gene_similarity_smooth],f)

"""
##This is where commentable section ends        
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
net_smooth_pickle = os.path.join(address_out_dir, 'Similarity_TCGA_RWR%s_%s' %(int(100 * restart_prob), network_name))
   
with open(net_smooth_pickle, 'rb') as f:
    nodes_genes_intersect, sample_name, expr_IQN_gene, gene_similarity_smooth  = pickle.load(f)  

###############################################################################
###############################################################################



####Focus only on genes shared between TCGA and GTEX
expr_IQN_GTEX = pd.read_csv(address_in_expr_GTEX, index_col=0)    #this is IQN of log2(RPKM+1) values for all cancers (not sorted)
pvalue_gt_df = pd.read_csv(address_out_pvalue_gt_generic, index_col=0, header=0) #Not sorted gene-TF correlation p-values
gene_TF_list_GTEX = list(expr_IQN_GTEX.index.values)

##############################################################################




for ci in range(len(cancer_list)):
    cancer = cancer_list[ci]
    pvalue_gp_array = np.zeros(len(nodes_genes_intersect))
    sample_interest = sample_name[ci]
    sample_other = []
    for cj in range(len(cancer_list)):
        if cj != ci:
            sample_other += sample_name[cj] #list of all samples that are not in cancer of interest
    for i in range(len(nodes_genes_intersect)):
        _, pvalue_gp_array[i] = ss.ttest_ind(expr_IQN_gene.loc[nodes_genes_intersect[i]][sample_interest], expr_IQN_gene.loc[nodes_genes_intersect[i]][sample_other], equal_var=False)
        if pvalue_gp_array[i] < eps:
            pvalue_gp_array[i] = eps

    #form the new p-values using Stouffer method:
    comb_gp_combined_list = [ss.combine_pvalues(pvalue_gp_array, method='stouffer', weights=gene_similarity_smooth[i,:]) for i in range(len(nodes_genes_intersect))] 

    ###Only focus on genes that are in both TCGA and GTEX:
    comb_gp_combined_dict = {nodes_genes_intersect[i]:comb_gp_combined_list[i] for i in range(len(nodes_genes_intersect)) if nodes_genes_intersect[i] in gene_TF_list_GTEX}
    gene_list_tmp = []
    stat_gp_combined_list_tmp = []
    pval_gp_combined_list_tmp = []
    for gene in comb_gp_combined_dict:
        gene_list_tmp.append(gene)
        stat_gp_combined_list_tmp.append(comb_gp_combined_dict[gene][0])
        if comb_gp_combined_dict[gene][1] < eps:
            pval_gp_combined_list_tmp.append(eps)
        else:   
            pval_gp_combined_list_tmp.append(comb_gp_combined_dict[gene][1])

    pval_gp_combined_array = np.array(pval_gp_combined_list_tmp)
    stat_gp_combined_array = np.array(stat_gp_combined_list_tmp)
    print(cancer, '# DEG %s' %(pval_gp_combined_array < 0.05).sum())

    #sort all genes based on pvalues and form dataframe:
    argsort_ind = np.argsort(stat_gp_combined_array)[::-1]
    pvalue_gp_array_psorted = pval_gp_combined_array[argsort_ind]
    only_gene_list_psorted = [gene_list_tmp[i] for i in argsort_ind]
    pvalue_gp_df_sorted = pd.DataFrame(pvalue_gp_array_psorted, index=only_gene_list_psorted, columns=['PValue'])   #this is the dataframe to be written to file

    #Now we form gene_TF dataframe where genes are sorted based on only_gene_list_psorted
    pvalue_gt_df_sorted = pvalue_gt_df.loc[only_gene_list_psorted]
    
    #write to file:
    address_pval_gp = os.path.join(address_out_dir, 'Pvalue_gene_phenotype_RWR%s_1vsAll_FPKM_%s.csv' %("{:.0f}".format(restart_prob*100), cancer))
    address_pval_tg = os.path.join(address_out_dir, 'Pvalue_TF_gene_RWR%s_%s_1vsAll_FPKM_%s.csv' %("{:.0f}".format(restart_prob*100), tf_gene_method, cancer))
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########
    pvalue_gp_df_sorted.to_csv(address_pval_gp)
    pvalue_gt_df_sorted.to_csv(address_pval_tg)
##########%%%%%%%%%%@@@@@@@@@@!!!!!!!!!!@@@@@@@@@@%%%%%%%%%%##########




