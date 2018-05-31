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

In order to run the code, you need to have Python 3.5 installed. In addition, the code uses the following python modules/libraries which need to be installed (the number in brackets show the version):
- [Numpy](http://www.numpy.org/) (version 1.13.0)
- [Scipy](https://www.scipy.org/) (verison 0.19.1)
- [Pandas](http://pandas.pydata.org/) (version 0.20.2)
- [Sklearn (scikit-learn)](http://scikit-learn.org/stable/) (version 0.18.1)
- [PyMC](https://pymc-devs.github.io/pymc/) (version 2.3.6)

Instead of installing all these libraries independently, you can use prebulit Python distributions such as [Anaconda](https://www.continuum.io/downloads), which provides a free academic subscription. If you are using Anaconda, you can easily install any specific version of the modules above using a command like:

conda install scipy=0.19.1

For pymc, you can use

conda install -c https://conda.binstar.org/pymc pymc=2.3.6

# Input files

### Description of required inputs:
#### Gene expression (features) file:
This is a genes x samples csv file where the first column contains name of genes and the first row contains name/IDs of the samples. ProGENI assumes that the expression of each gene (across all samples) follows a normal distribution. As a result, we recommend you perform proper transformation on your expression data (e.g. log2 transform on microarray data) to satsify this condition for best results. NAs are not allowed in this file. 

Example Gene expression file:

|  | sample_1 | sample_2 | sample_3 |
| :--- | :--- | :--- | :--- |
| G1 | 0.24 | 0.67 | 2.12 |  
| G2 | 0.34 | -1.34 | 0.45 |
| G3 | 1.51 | 0.05 | -0.22 |
| G4 | 0.03 | 0.55 | 1.15 |
| G5 | -0.23 | 0.23 | 0.55 |
| G6 | 0.94 | 0.33 | 1.12 |


#### Phenotype (response) file:
This is a phenotype x samples csv file where the first column contains name of different phenotypes (e.g. different drugs) and the first row contains name/IDs of the samples. Make sure that the samples are ordered the same as the gene expression file (i.e. the first row of both files hould be identical). NAs are allowed in this file. 

Example phenotype file:

|  | sample_1 | sample_2 | sample_3 |
| :--- | :--- | :--- | :--- |
| drug_1 | 0.65 | 0.12 | 1.45 |  
| drug_2 | 1.67 | NA | 2.45 |
| drug_3 | 2.51 | 0.56 | 0.34 |


#### Network edge file:
This is a csv file which contains information on gene-gene interactions. The first should be the header of the file. The network should be represented as a three-column format where each edge in the network is represented as a row in the file: the first two columns contain name of genes and the third column shows the (positive) weight (e.g. representing confidence) corresponding to this relationship. If the set of genes in the network is slightly different from the set of genes in the gene expression data, ProGENI will focus on the intersection of the genes.  

Example network edge file:

| node_1 | node_2 | weight |
| :--- | :--- | :--- |
| G1 | G4 | 777 |
| G1 | G6 | 232 |
| G2 | G4 | 999 |
| G2 | G7 | 131 |
| G4 | G5 | 444 |
 

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
