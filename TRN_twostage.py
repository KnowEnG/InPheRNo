# -*- coding: utf-8 -*-
"""
This script uses the twostage model to learn TRN on real data. 
"""

import numpy as np
from pymc import Beta, Uniform, DiscreteUniform, deterministic
from pymc import Exponential, Poisson, Uniform, Bernoulli
import pymc
import os
import pickle
import argparse
import pandas as pd
import sys
from scipy.stats import beta, uniform
from scipy.stats import rv_continuous

###############################################################################
def parse_args():
    """
    Parse the arguments.
    Parse the command line arguments/options using the argparse module
    and return the parsed arguments (as an argparse.Namespace object,
    as returned by argparse.parse_args()).
    Returns:
        argparse.Namespace: the parsed arguments
    """
    parser = argparse.ArgumentParser()
#    parser.add_argument('-l', '--location', default='cloud9', help='laptop or cloud9')
#    parser.add_argument('-s', '--seed', default=1011, help='seed used for random generator')
    parser.add_argument('-atg0', '--A_TF_gene_h0', default="None", help='alpha in beta dist. of TF-gene under Null hypothesis (TF does not regulate gene). If None, it is learnt from data. Otherwise a number must be given.')
    parser.add_argument('-atg1', '--A_TF_gene_h1', default="None", help='alpha in beta dist. of TF-gene under alternative hypothesis (TF regulates gene). If None, it is learnt from data. Otherwise a number must be given.')
    parser.add_argument('-agp', '--A_gene_pheno', help='alpha in beta dist. of gene-pheno. It is a number and must be given.')
    parser.add_argument('-pt', '--Prior_T', default="None", help='Prior probability of T=1')
    parser.add_argument('-ptm', '--Prior_T_method', default='fixed', help='Method to set prior probability of T=1. Options are fixed and learn. If fixed is chosen and args.priot_T is None, we use 1/(10*n_TF).')
    parser.add_argument('-rtg', '--R_TF_gene', default="None", help='mixing parameter')
    parser.add_argument('-igp', '--input_gene_pheno', default='./Results/Pvalue_gene_phenotype_interest_tmp.csv', help='Full address of file containing gene-phenotype association p-values for genes of interest')
    parser.add_argument('-itg', '--input_TF_gene', default='./Results/Pvalue_gene_tf_tmp.csv', help='Full address of file containing TF-gene correlation p-values')
    parser.add_argument('-od', '--output_dir', default='./tmp', help='address of directory for results')
    parser.add_argument('-of', '--output_file', default="None", help='prefix to be added to output files')
    parser.add_argument('-mnt', '--max_num_TF', default=15, help='Max number of TFs identified using first stage')
    parser.add_argument('-ir', '--index_repeat', default=100, help='Index of the repeat. Repeats are used to ensure stability of results.')
    parser.add_argument('-ni', '--num_iteration', default=200, help='number of iterations for the PGM')
    parser.add_argument('-nb', '--num_burn', default=100, help='number of iterations to burn')
    parser.add_argument('-nt', '--num_thin', default=1, help='number of iterations to thin by')
    parser.add_argument('-si', '--start_index', default=0, help='start index of the gene to consider')
    parser.add_argument('-ei', '--end_index', default='None', help='end index of the gene to consider')

    args = parser.parse_args()
    return(args)
###############################################################################

def Model_twostage_fit_v2(n_TF, n_gene, p_gene_array, p_TF_gene_array, num_iter, num_burn, num_thin, prior_T, prior_T_method, r_TF_gene, a_TF_gene_h1, a_TF_gene_h0, a_gene):
    """
    Assumptions: We allow learning of parameters
    """
    a_gp = float(a_gene)
    if a_TF_gene_h0 == 'None':
        a_tg0 = Uniform('a_tg0', lower=0.5, upper=1)
    else:
        a_tg0 = float(a_TF_gene_h0)
    if a_TF_gene_h1 == 'None':
        a_tg1 = Uniform('a_tg1', lower=0, upper=0.5)
    else:
        a_tg1 = float(a_TF_gene_h1)
    p_T = float(prior_T)

    if r_TF_gene == 'None':
        r_tg = Uniform('r_tg', lower=0, upper=1)
    else:
        r_tg = float(r_TF_gene)
    p_gene = np.zeros(n_gene, dtype=object)     #the ovserved variables 
    T = np.zeros((n_TF, n_gene), dtype=object)  #variables showing TF-gene-pheno relationship
    T_sum = np.zeros(n_gene, dtype=object)
    p_TF_gene = np.zeros((n_TF, n_gene), dtype=object)  #p-value of correlation of gene TF
    for j in range(n_gene):
        for i in range(n_TF):
            T[i, j] = Bernoulli('T_%i_%i' %(i, j), p=p_T)
            
            #If T[i, j] = 0: then p_TF_gene is coming from a mixture of beta and uniform (r is the mixture param)
            @pymc.stochastic(name='p_TF_gene_%i_%i' %(i, j), dtype=float, observed=True)
            def temp_p_TF_gene(value=p_TF_gene_array[i, j], TF_gene_ind=T[i, j], a0=a_tg0, a1=a_tg1, r=r_tg) :
                if TF_gene_ind:
                    out = pymc.distributions.beta_like(value, alpha=a1, beta=1)
                else:
                    out = np.log(r * np.exp(pymc.distributions.beta_like(value, alpha=a1, beta=1))
                                + (1 - r) * np.exp(pymc.distributions.beta_like(value, alpha=a0, beta=1)))
                return out
            p_TF_gene[i, j] = temp_p_TF_gene
            
        #we define a deterministic function to find values of T           
        @pymc.deterministic(name='T_sum_%i' %j, plot=False)
        def temp_T_sum(ind_vec=T[:,j]): 
            return (np.sum(ind_vec)>0)
        T_sum[j] = temp_T_sum

        #If T_sum[j] == 0: then p_TF_gene is coming from a uniform; else, beta
        @pymc.stochastic(name='p_gene_%i' %j, dtype=float, observed=True)
        def temp_p_gene(value=p_gene_array[j], ind=T_sum[j], a=a_gp):
            if ind:
                out = pymc.distributions.beta_like(value, alpha=a, beta=1)
            else:
                out = pymc.distributions.uniform_like(value, 0, 1)
            return out
        p_gene[j] = temp_p_gene
    if a_gene == None and a_TF_gene_h0 == None and a_TF_gene_h1 == None:
        M5 = pymc.MCMC([T, T_sum, a_gp, a_tg0, a_tg1])
    else:        
        M5 = pymc.MCMC([T, T_sum])
    M5.sample(iter=int(num_iter), burn=int(num_burn), thin=int(num_thin))
    return(M5)



###############################################################################
args = parse_args()
#address_gene_pheno = os.path.join(args.input_dir, args.input_gene_pheno)
address_gene_pheno = args.input_gene_pheno
address_TF_gene = args.input_TF_gene



delim_gp = ','
delim_tg = ','
if args.input_gene_pheno[-3:] in ['tsv', 'txt']:
    delim_gp = '\t'
if args.input_TF_gene[-3:] in ['tsv', 'txt']:
    delim_tg = '\t'

pvalue_gene_pheno = pd.read_csv(address_gene_pheno, sep=delim_gp, index_col=0, header=0).T
pvalue_TF_gene = pd.read_csv(address_TF_gene, sep=delim_gp, index_col=0, header=0).T

if (pvalue_gene_pheno.columns != pvalue_TF_gene.columns).sum() > 0: 
    sys.exit('Input files do not match!')
else:
    print('Matching inputs!')


n_TF = len(pvalue_TF_gene.index)
n_gene = len(pvalue_TF_gene.columns)

eps = 1e-150 #we cap the minimum p-values to eps to avoid precision issues.
pvalue_TF_gene[(pvalue_TF_gene<eps) & (pvalue_TF_gene>-0.5)] = eps
pvalue_gene_pheno[pvalue_gene_pheno < eps] = eps




################################################################################
A_TF_gene_h1 = args.A_TF_gene_h1  #beta distribution parameter of alternative hypothesis. If None, it will be learnt.
A_TF_gene_h0 = args.A_TF_gene_h0  #beta distribution parameter of Null hypothesis. If None, it will be learnt.
A_gene = args.A_gene_pheno     #beta distribution parameter for gene phenotype association. 


print('A_gene_pheno = ', A_gene)

R_tg = args.R_TF_gene

Num_iter = args.num_iteration
Num_burn = args.num_burn
Num_thin = args.num_thin
Prior_T = args.Prior_T 
if Prior_T == "None":
    Prior_T = 1/int(args.max_num_TF) * 0.1
    print('Prior_T = ', Prior_T)

R_TF_gene = args.R_TF_gene
Prior_T_method = args.Prior_T_method
#np.random.seed(int(args.seed))


end_index = args.end_index
if end_index == "None":
    end_index = n_gene
    

if args.output_file == "None":
    address_outputfile = os.path.join(args.output_dir, 'InPheRNo_tmp_out_repeat%s.csv' %(args.index_repeat))        
else:
    address_outputfile = os.path.join(args.output_dir, args.output_file + 'InPheRNo_tmp_out_repeat%s.csv' %(args.index_repeat))        

if int(args.start_index) == 0:
    T_recovered = np.zeros((n_TF, n_gene))   #A matrix that shows probabilities of each variable being 1 
    print(np.shape(pvalue_TF_gene))
else:
    T_recovered_DF = pd.read_csv(address_outputfile, index_col=0, header=0).T
    T_recovered = T_recovered_DF.values

################################################################################

for jj in range(int(args.start_index), int(end_index)): 
    print('\n', 'repeat', args.index_repeat, 'gene', jj)
    selected_ind = pvalue_TF_gene.iloc[:, jj].values > -0.5
    n_selected_TF = len(pvalue_TF_gene.iloc[selected_ind, jj][:, None])
    model = Model_twostage_fit_v2(n_selected_TF, 1, np.asarray([pvalue_gene_pheno.iloc[0][jj]]), pvalue_TF_gene.iloc[selected_ind, jj][:, None], Num_iter, Num_burn, Num_thin, Prior_T, Prior_T_method, R_TF_gene, A_TF_gene_h1, A_TF_gene_h0, A_gene)
    T_recovered_tmp = np.zeros((n_selected_TF))
    num_samp = len(model.trace('T_0_0')[:])
    for i in range(n_selected_TF):
        T_recovered_tmp[i] = model.trace('T_%i_%i' %(i, 0))[:].sum() / num_samp
    T_recovered[selected_ind, jj] = T_recovered_tmp
    #print(pvalue_TF_gene.iloc[selected_ind, jj][:, None])
    #print(T_recovered_tmp)

T_recovered_out_DF = pd.DataFrame(T_recovered.T, index=pvalue_TF_gene.columns, columns=pvalue_TF_gene.index)
T_recovered_out_DF.to_csv(address_outputfile)