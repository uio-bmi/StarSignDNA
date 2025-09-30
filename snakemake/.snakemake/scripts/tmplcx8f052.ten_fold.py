######## snakemake preamble start (automatically inserted, do not edit) ########
import sys;sys.path.extend(['/home/chbope/miniconda3/lib/python3.13/site-packages', '/home/chbope/Documents/starsignDNA/StarSignDNA_srcode/StarSignDNA_srcode/StarSignDNA/snakemake', '/home/chbope/miniconda3/bin', '/home/chbope/miniconda3/lib/python3.13', '/home/chbope/miniconda3/lib/python3.13/lib-dynload', '/home/chbope/miniconda3/lib/python3.13/site-packages', '/home/chbope/Documents/starsignDNA/StarSignDNA_srcode/StarSignDNA_srcode/StarSignDNA', '/home/chbope/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpqjzf7uav/file/home/chbope/Documents/starsignDNA/StarSignDNA_srcode/StarSignDNA_srcode/StarSignDNA/snakemake/script', '/home/chbope/Documents/starsignDNA/StarSignDNA_srcode/StarSignDNA_srcode/StarSignDNA/snakemake/script']);import pickle;from snakemake import script;script.snakemake = pickle.loads(b'\x80\x04\x95\x12\x04\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c\x1aresults/data/4/train_m.csv\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x10h\x06\x8c\x0eAttributeGuard\x94\x93\x94)\x81\x94}\x94\x8c\x04name\x94h\x10sbh\x11h\x13)\x81\x94}\x94h\x16h\x11sbub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x1bresults/data/4/11/0.2/e.csv\x94\x8c\x1bresults/data/4/11/0.2/s.csv\x94e}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x13)\x81\x94}\x94h\x16h\x10sbh\x11h\x13)\x81\x94}\x94h\x16h\x11sbub\x8c\r_params_store\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x13)\x81\x94}\x94h\x16h\x10sbh\x11h\x13)\x81\x94}\x94h\x16h\x11sbub\x8c\r_params_types\x94}\x94\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x04data\x94\x8c\x014\x94\x8c\x0211\x94\x8c\x030.2\x94e}\x94(h\x0c}\x94(\x8c\x07dataset\x94K\x00N\x86\x94\x8c\x04fold\x94K\x01N\x86\x94\x8c\x01k\x94K\x02N\x86\x94\x8c\x06Lambda\x94K\x03N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x13)\x81\x94}\x94h\x16h\x10sbh\x11h\x13)\x81\x94}\x94h\x16h\x11sb\x8c\x07dataset\x94h7\x8c\x04fold\x94h8hAh9\x8c\x06Lambda\x94h:ub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x13)\x81\x94}\x94h\x16h\x10sbh\x11h\x13)\x81\x94}\x94h\x16h\x11sbhUK\x01hWK\x01hYhRub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x13)\x81\x94}\x94h\x16h\x10sbh\x11h\x13)\x81\x94}\x94h\x16h\x11sbub\x8c\x06config\x94}\x94\x8c\x04rule\x94\x8c\x02cv\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8ce/home/chbope/Documents/starsignDNA/StarSignDNA_srcode/StarSignDNA_srcode/StarSignDNA/snakemake/script\x94ub.');del script;from snakemake.logging import logger;from snakemake.script import snakemake;__real_file__ = __file__; __file__ = '/home/chbope/Documents/starsignDNA/StarSignDNA_srcode/StarSignDNA_srcode/StarSignDNA/snakemake/script/ten_fold.py';
######## snakemake preamble end #########
#!/usr/bin/env python3
"""
Ten-fold cross-validation script for mutational signature analysis.

This script performs signature extraction on training data using specified
hyperparameters (number of signatures k and regularization lambda).
"""

import pandas as pd
import numpy as np
from mutation_sim_train import main

# Extract hyperparameters from snakemake wildcards
k = int(snakemake.wildcards.k)  # Number of signatures
Lambda = float(snakemake.wildcards.Lambda)  # Regularization parameter

# Read the training mutation matrix (exclude the first column which is likely sample IDs)
train_m = pd.read_csv(snakemake.input[0], delimiter="\t")
train_m = train_m.iloc[:, 1:]  # Remove first column (sample IDs)

# Convert to numpy array for processing
train_m = train_m.to_numpy().astype(float)

# Run signature extraction to get exposures and signatures
E, S = main(train_m, k, Lambda)

# Save results
np.savetxt(snakemake.output[0], E)  # Save exposures
np.savetxt(snakemake.output[1], S)  # Save signatures
