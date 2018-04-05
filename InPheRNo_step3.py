
"""
This code combines the posterior probability obtained using PGM from different
repeats.
"""

import pandas as pd
import numpy as np
import os
import argparse



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
    parser.add_argument('-id', '--input_dir', default='./tmp', help='address of directory for results of different repeats')
    parser.add_argument('-if', '--input_file', default="None", help='prefix of the repeat files')
    parser.add_argument('-nr', '--num_repeat', default=100, help='Number of the repeats. Repeats are used to ensure stability of results. At least 100 repeats are recommended.')
    parser.add_argument('-od', '--output_dir', default='./Results', help='address of output directory')
    parser.add_argument('-on', '--output_network', default = 'Final_phenotype_relevant_TRN.csv', help = 'The final result. The scores assigned to each edge can be thresholded to obtain a network (e.g. >0.5).')

    args = parser.parse_args()
    return(args)


###############################################################################
args = parse_args()

num_repeat = int(args.num_repeat) #number of repetitions to combine 


address_mean_out = os.path.join(args.output_dir, args.output_network)

###############################################################################
##read repeat 0 for initialization:
if args.input_file in ["None", "NONE"]:
    address_in0 = os.path.join(args.input_dir, 'InPheRNo_tmp_out_repeat0.csv')        
else:
    address_in0 = os.path.join(args.input_dir, args.input_file + 'InPheRNo_tmp_out_repeat0.csv')        

posterior0 = pd.read_csv(address_in0, index_col=0, header=0)

all_posterior = np.zeros((num_repeat, np.shape(posterior0)[0], np.shape(posterior0)[1])) # a numpy array where axis 0 shows repeats


for i in range(10): #range(num_repeat):
    if args.input_file == "None":
        address_inputfile = os.path.join(args.input_dir, 'InPheRNo_tmp_out_repeat%s.csv' %i)        
    else:
        address_inputfile = os.path.join(args.input_dir, args.output_file + 'InPheRNo_tmp_out_repeat%s.csv' %i)        

    postrior_in = pd.read_csv(address_inputfile, index_col=0, header=0)
    all_posterior[i, :, :] = postrior_in.values


all_posterior_mean = np.mean(all_posterior, axis=0)
all_posterior_mean /= np.max(all_posterior_mean)

all_mean_DF = pd.DataFrame(all_posterior_mean, index=posterior0.index, columns=posterior0.columns)

all_mean_DF.to_csv(address_mean_out)

