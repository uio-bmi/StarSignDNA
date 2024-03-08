#!/usr/bin/env python3
import pandas as pd
from mutation_sim_test import main
import numpy as np
snakemake = snakemake
test_m = pd.read_csv(snakemake.input[0], delimiter="\t")
test_m = test_m.iloc[:,1:]
###
#headers = test_m.columns
#mutation_counts = test_m.to_dict(orient='list')
#mutation_counts_df = pd.DataFrame(mutation_counts)
#test_m  = get_tri_context_fraction(mutation_counts_df)
#print("TTTTT",test_m)
###
test_m = test_m.to_numpy().astype(float)
S = pd.read_csv(snakemake.input[1], header=None, delimiter=" ")
E, loss = main(test_m, S)

np.savetxt(snakemake.output[0], E)
open(snakemake.output[1], "w").write(str(loss))
