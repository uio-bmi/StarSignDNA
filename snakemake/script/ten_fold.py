#!/usr/bin/env python3
import pandas as pd
import numpy as np
from mutation_sim_train import main
snakemake = snakemake

k = int(snakemake.wildcards.k)
Lambda = float(snakemake.wildcards.Lambda)

train_m = pd.read_csv(snakemake.input[0], delimiter="\t")
train_m = train_m.iloc[:,1:]
###
#headers = train_m.columns
#mutation_counts = train_m.to_dict(orient='list')
# Convert mutation counts to a DataFrame
#mutation_counts_df = pd.DataFrame(mutation_counts)
#print("Mut", mutation_counts_df)
#print(M)
#M = M.to_numpy().astype(float)
#train_m  = get_tri_context_fraction(mutation_counts_df)
###
train_m = train_m.to_numpy().astype(float)
E, S = main(train_m, k, Lambda)

np.savetxt(snakemake.output[0], E)
np.savetxt(snakemake.output[1], S)
