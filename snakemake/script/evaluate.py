#!/usr/bin/env python3
import pandas as pd
from .refit import refit
import numpy as np
snakemake = snakemake
test_m = pd.read_csv(snakemake.input[0], header=None, delimiter="\t")
S = pd.read_csv(snakemake.input[1], header=None, delimiter=" ")
E, loss = refit(test_m, S)

np.savetxt(snakemake.output[0], E)
open(snakemake.output[1], "w").write(str(loss))
