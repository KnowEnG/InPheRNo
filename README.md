# *InPheRNo* - Inference of Phenotype-relevant Regulatory Networks
#### Amin Emad (email: emad2 (at) illinois (dot) edu)
#### KnowEnG BD2K Center of Excellence
#### University of Illinois Urbana-Champaign


# Motivation
Reconstruction of transcriptional regulatory networks (TRNs) is a powerful approach to unravel the gene expression programs involved in healthy and disease states of a cell. However, these networks are usually reconstructed independent of the phenotypic properties of the samples and therefore cannot identify regulatory mechanisms that are related to a phenotypic outcome of interest. InPheRNo (Inference of Phenotype-relevant Regulatory Networks) is a computational tool to reconstruct phenotype-relevant transcriptional regulatory networks (TRNs) using transcriptomic data.
Here we present InPheRNo, a novel computational tool to reconstruct ‘phenotype-relevant’ transcriptional regulatory networks. This method is based on a probabilistic graphical model whose conditional probability distributions model the simultaneous effects of multiple transcription factors (TFs) on their target genes as well as the statistical relationship between target gene expression and phenotype. 


The figure below depcits the method overview. 

![Method Overview](/Images/Figure1_Pipeline_final.png)

# Requirements

In order to run the code, you need to have Python 3.5 installed. In addition, the code uses the following python modules/libraries which need to be installed (the number in brackets show the version of the module used to generate the results in the manuscript):
- [Numpy](http://www.numpy.org/) (version 1.13.0)
- [Scipy](https://www.scipy.org/) (verison 0.19.1)
- [Pandas](http://pandas.pydata.org/) (version 0.20.2)
- [Sklearn (scikit-learn)](http://scikit-learn.org/stable/) (version 0.18.1)
- [PyMC](https://pymc-devs.github.io/pymc/) (version 2.3.6)

Instead of installing all these libraries independently, you can use prebulit Python distributions such as [Anaconda](https://www.continuum.io/downloads), which provides a free academic subscription. If you are using Anaconda, you can easily install any specific version of the modules above using a command like:

conda install pymc=2.3.6

# Running InPheRNo
Running InPheRNo involves running three manuscripts in a row. Since the intermediate results are used in the following steps, one needs to wait for the preceeding step to finish before running the next step. 

## STEP 1:
### Description of required inputs:
#### Input1: A file containing the list of transcription factors (TFs):
This is a csv file in which rows contain the names of the regulators (e.g. TFs). The file should not have a header. As an example see the file "Data/TF_Ensemble.csv". 

#### Input2: A file containing p-values of gene-phenotype associations only for genes of interest:
This is a (gene x phenotype) csv file (see "Data/Pvalue_gene_phenotype_interest.csv" as an example). The rows correspond to target genes of interest (this may be only a subset of all genes, or it may be all the genes). The p-value for TF-phenotype should not be included in this file. The value assigned to each gene represents the p-value of association between the expression of that gene and the variation in the phenotype across different samples obtained using a proper statistical test (e.g. a ttest for binary phenotype or Pearson's correlation for continuous, etc.). The genes should be sorted in an ascending order based on the p-value (smallest p-values appear first). The file is assumed to have a header. 

Example:

|  | Pvalue |
| :--- | :--- |
| gene1 | 1E-22 |  
| gene2 | 5E-14 |
| gene3 | 3E-10 |


#### Input3: A file containing gene and TF expression data.'):
This is a (gene x samples) csv file containing the normalized gene (and TF) expression profiles across different samples. This file must contain expression of target genes provided in Input2 and TFs provided in Input1. The file has a header representing sample names. See "Data/expr_sample.csv" as a sample input.  

Example:

| sample1 | sample2 | sample3 |
| :--- | :--- | :--- |
| TF1 | 0.1 | 0.9 |
| TF2 | -0.3 | 0.5 |
| gene1 | -1.1 | 0.6 |
| gene2 | 0.9 | -2.3 |
| gene3 | 0.4 | 0.8 |
 

# Running ProGENI
### With default settings
There are only 3 required (positional) arguments that needs to be specified by the user:
- input_expression: name of the csv file containing the gene expression data
- input_response: name of the csv file containing the phenotype data
- input_network: name of the csv file containing the network edges
By default, ProGENI assumes that all these files are located in the current directory. Given these arguments, one can run ProGENI with default settings. The results will be saved in a file called "results.csv" in the current directory, which contains the ranked list of genes for each response. Only genes shared between the network and gene expression data will be included in the results. The following line shows how to run ProGENI:
```
python3 ProGENI.py gene_expr.csv phenotype.csv network.csv
```

### With advanced settings
In addition to the positional arguemtns, one can use the following optional arguments to change the default settings.
- -o, --output (string, default='results.csv'): name of the file containg the results
- -do, --directory_out (string, default='./'): directory for the results
- -de, --directory_expression (string, default='./'): directory containing the gene expression file
- -dr, --directory_response (string, default='./'): directory containing the response file
- -dn --directory_network (string, default='./'): directory containing the network file
- -s, --seed (integer, default=1011): the seed for the pseudo random generator used in bootstrap sampling
- -nr, --num_RCG (integer, default=100): number of genes in the response-correlated gene (RCG) set
- -pt, --prob_restart_trans (float, default=0.5): restart probability of RWR to network-transform gene expression
- -pr, --prob_restart_rank (float, default=0.5): restart probability for RWR used to rank nodes w.r.t. RCG
- -t, --tolerance (float, default=1e-8): residual tolerance used to determine convergence of RWR
- -mi, --max_iteration (integer, default=100): maximum number of iterations used in RWR
- -nb, --num_bootstrap (integer, default=1): number of bootstrap samplings
- -pb, --percent_bootstrap (integer, default=100): percent of samples for bootstrap sampling (between 0-100)

For example, to run Robust-ProGENI with 80% bootstrap sampling and 50 times repeat and save the results in a file called "results_80_50.csv" one can use the following line:
```
python3 ProGENI.py gene_expr.csv phenotype.csv network.csv -o results_80_50.csv -nb 50 -pb 80
```

### Sample inputs and outputs:
To test whether ProGENI runs as expected on your machine, you can use the sample inputs in the folder "sample_data" and run ProGENI with num_RCG=2 and other arguments set with default values. The results should match the file "results_sample.csv".
```
python3 ProGENI.py gene_expr_sample.csv response_sample.csv network_sample.csv -nr 2
```
# Running ProGENI_simplified
ProGENI_simplified.py provides a simplified implementation of ProGENI. In this variation, the Pearson correlation coefficient of network transformed gene expressions and phenotype is used to rank the genes. In other words, the steps involving identification of a RCG set and ranking genes in the network with respect to the RCG are removed. This method is called "ProGENI-PCC" in the manuscript. 
Usage of this variation is very similar to ProGENI.py, except that -nr and -pr do not need to be provided. 
