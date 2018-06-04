# *InPheRNo* - Inference of Phenotype-relevant Regulatory Networks
#### Amin Emad (email: emad2 (at) illinois (dot) edu)
#### KnowEnG BD2K Center of Excellence
#### University of Illinois Urbana-Champaign


# Motivation
Reconstruction of transcriptional regulatory networks (TRNs) is a powerful approach to unravel the gene expression programs involved in healthy and disease states of a cell. However, these networks are usually reconstructed independent of the phenotypic properties of the samples and therefore cannot identify regulatory mechanisms that are related to a phenotypic outcome of interest. InPheRNo (Inference of Phenotype-relevant Regulatory Networks) is a computational tool to reconstruct phenotype-relevant transcriptional regulatory networks (TRNs) using transcriptomic data.
Here we present InPheRNo, a novel computational tool to reconstruct ‘phenotype-relevant’ transcriptional regulatory networks. This method is based on a probabilistic graphical model (PGM) whose conditional probability distributions model the simultaneous effects of multiple transcription factors (TFs) on their target genes as well as the statistical relationship between target gene expression and phenotype. 


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
Running InPheRNo involves running three manuscripts (InPheRNo_step1.py, InPheRNo_step2.py and InPheRNo_step3.py) in a row. Since the intermediate results are used in the following steps, one needs to wait for the preceeding step to finish before running the next step. 

## STEP 1:
### Description of required inputs:
#### Input1.1: A file containing the list of transcription factors (TFs):
This is a csv file in which rows contain the names of the regulators (e.g. TFs). The file should not have a header. As an example see the file "Data/TF_Ensemble.csv". 

#### Input1.2: A file containing p-values of gene-phenotype associations only for genes of interest:
This is a (gene x phenotype) csv file (see "Data/Pvalue_gene_phenotype_interest.csv" as an example). The rows correspond to target genes of interest (this may be only a subset of all genes, or it may be all the genes). The p-value for TF-phenotype should not be included in this file. The value assigned to each gene represents the p-value of association between the expression of that gene and the variation in the phenotype across different samples obtained using a proper statistical test (e.g. a ttest for binary phenotype or Pearson's correlation for continuous, etc.). The genes should be sorted in an ascending order based on the p-value (smallest p-values appear first). The file is assumed to have a header. 

Example:

|  | Pvalue |
| :--- | :--- |
| gene1 | 1E-22 |  
| gene2 | 5E-14 |
| gene3 | 3E-10 |

#### Input1.3: A file containing gene and TF expression data:
This is a (gene x samples) csv file containing the normalized gene (and TF) expression profiles across different samples. This file must contain expression of target genes provided in Input1.2 and TFs provided in Input1.1. The file has a header representing sample names. See "Data/expr_sample.csv" as a sample input.  

Example:

|  | sample1 | sample2 | sample3 |
| :--- | :--- | :--- | :--- |
| TF1 | 0.1 | 0.9 | 0.5 |
| TF2 | -0.3 | 0.5 | -0.6 |
| gene1 | -1.1 | 0.6 | 1.4 |
| gene2 | 0.9 | -2.3 | -0.3 |
| gene3 | 0.4 | 0.8 | 1.5 |
 
### Description of outputs:
The first step of InPheRNo generates two output files that by default will be located in a directory called "Results" placed in the current directory. These intermediate outupts will be used in the next step of InPheRNo. 

#### Output1.1: Pvalue_gene_phenotype_interest_tmp.csv
If default parameters are used to run the first step, Output1.1 will be a file called "Pvalue_gene_phenotype_interest_tmp.csv" which is generated from Input1.2, properly sorted and cleaned up (if necessary). See folder "Results" for a sample.

#### Output1.2: Pvalue_gene_tf_tmp.csv
If default parameters are used to run the first step, Output1.2 will be a file called "Pvalue_gene_tf_tmp.csv". This is a (gene x TF)  csv file containing p-values of gene-tf association, sorted in an ascending order based on Output1.1 file. The file has a header. See folder "Results" for a sample.

### Running InPheRNo_step1.py: 
#### With default settings
To Run this step with default parameters, place all the three input files above in one folder. Then specify the following four arguments:
- input_directory: address of the data directory containing the three input files (e.g. "./Data")
- input_tf: name of Input1.1 file containing the name of regulators (e.g. "TF_Ensemble.csv")
- input_gene_phenotype_interest: name of Input1.2 containing p-value of gene-phenotype (e.g. "Pvalue_gene_phenotype_interest.csv")
- input_expression: name of Input1.3 containing the expression of genes and TFs (e.g. "expr_sample.csv")

The following line shows how to run InPheRNo using the sample files:
```
python3 InPheRNo_step1.py --input_directory ./Data --input_tf TF_Ensemble.csv --input_gene_phenotype_interest Pvalue_gene_phenotype_interest.csv --input_expression expr_sample.csv
```

By default, InPheRNo writes the intermediate outputs generated in this step into a directory called "Results" in the current directory. To change the location of the intermediate results, see advanced settings. 

#### With advanced settings
In addition to the arguments above, one can use the following optional arguments to change the default settings.
- -od, --output_directory (string, default='./Results'): Address of the output directory
- -mt, --max_num_tf (integer, default = 15): Maximum number of TFs recovered for each gene using Elastic Net
- -lr, --l1_ratio, (float, default = 0.5): l1 ratio of the Elastic Net model
- -ogp, --output_gene_phenotype (string, default = 'Pvalue_gene_phenotype_interest_tmp.csv'): Name of Output1.1 file
- -tgt, --output_gene_tf (string, default = 'Pvalue_gene_tf_tmp.csv'): Name of Output1.2 file


## STEP 2:
### Description of the required inputs:
This step requires three input files. Two of these input files are the intermediate outputs generated in STEP1.

#### Input2.1: A file containing p-values of gene-phenotype associations for all the genes:
This is a (gene x phenotype) csv file (see "Data/Pvalue_gene_phenotype_all.csv" as an example) and is very similar to Input 1.2. The main difference is that this file needs to contain the gene-phenotype association p-values for all the genes and not just the genes of interest. This file is used to estimate the parameters of the PGM. The rows correspond to target genes of interest (this may be only a subset of all genes, or it may be all the genes). The p-value for TF-phenotype should not be included in this file. The value assigned to each gene represents the p-value of association between the expression of that gene and the variation in the phenotype across different samples obtained using a proper statistical test (e.g. a ttest for binary phenotype or Pearson's correlation for continuous, etc.). The genes should be sorted in an ascending order based on the p-value (smallest p-values appear first). The file is assumed to have a header. 

#### Input2.2: The intermediate Output1.1
This is the Output1.1 file generated in STEP1. If default parameters are used to run the first step, this will be a file called "Pvalue_gene_phenotype_interest_tmp.csv" which is generated from Input1.2, properly sorted and cleaned up (if necessary). See folder "Results" for a sample.

#### Input2.3: The intermediate Output1.2
This is the Output1.2 file generated in STEP1. If default parameters are used to run the first step, Output1.2 will be a file called "Pvalue_gene_tf_tmp.csv". This is a (gene x TF) csv file containing p-values of gene-tf association, sorted in an ascending order based on Output1.1 file. The file has a header. See folder "Results" for a sample.












### Sample inputs and outputs:
To test whether ProGENI runs as expected on your machine, you can use the sample inputs in the folder "sample_data" and run ProGENI with num_RCG=2 and other arguments set with default values. The results should match the file "results_sample.csv".
```
python3 ProGENI.py gene_expr_sample.csv response_sample.csv network_sample.csv -nr 2
```
# Running InPheRNo_simplified
ProGENI_simplified.py provides a simplified implementation of ProGENI. In this variation, the Pearson correlation coefficient of network transformed gene expressions and phenotype is used to rank the genes. In other words, the steps involving identification of a RCG set and ranking genes in the network with respect to the RCG are removed. This method is called "ProGENI-PCC" in the manuscript. 
Usage of this variation is very similar to ProGENI.py, except that -nr and -pr do not need to be provided. 
