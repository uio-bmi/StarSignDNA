#!/usr/bin/env python3
import pandas as pd
import numpy as np
from cucumber import denovo
snakemake = snakemake

k = int(snakemake.wildcards.k)
Lambda = float(snakemake.wildcards.Lambda)

train_m = pd.read_csv(snakemake.input[0], header=None, delimiter="\t")
train_m = train_m.iloc[:,1:]
E, S = denovo(train_m, k, Lambda)

np.savetxt(snakemake.output[0], E)
np.savetxt(snakemake.output[1], S)
