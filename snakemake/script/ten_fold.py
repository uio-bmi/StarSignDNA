#!/usr/bin/env python3
import pandas as pd
import numpy as np
from mutation_sim_train import main
snakemake = snakemake

k = int(snakemake.wildcards.k)
Lambda = float(snakemake.wildcards.Lambda)

train_m = pd.read_csv(snakemake.input[0],header=None, delimiter="\t")
train_m = train_m.iloc[:,1:]
train_m = train_m.to_numpy().astype(float)
E, S = main(train_m, k, Lambda)

np.savetxt(snakemake.output[0], E)
np.savetxt(snakemake.output[1], S)
