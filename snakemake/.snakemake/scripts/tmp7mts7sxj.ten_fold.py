
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/opt/homebrew/lib/python3.9/site-packages', '/Users/bope/Documents/MutSig/scientafellow/packages/Projet/mutation_L1/python/m_results/manuscript/package/cucumber/snakemake/script']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xd8\x04\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c\x1aresults/data/3/train_m.csv\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x10\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x16)}\x94\x8c\x05_name\x94h\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x19results/data/3/10/0/e.csv\x94\x8c\x19results/data/3/10/0/s.csv\x94e}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x04data\x94\x8c\x013\x94\x8c\x0210\x94\x8c\x010\x94e}\x94(h\x0c}\x94(\x8c\x07dataset\x94K\x00N\x86\x94\x8c\x04fold\x94K\x01N\x86\x94\x8c\x01k\x94K\x02N\x86\x94\x8c\x06Lambda\x94K\x03N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94b\x8c\x07dataset\x94hD\x8c\x04fold\x94hEhNhF\x8c\x06Lambda\x94hGub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c0/var/folders/21/n1m80vz54z3djk8bv2d22z940000gn/T\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bhfK\x01hhK\x01hjhcub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06config\x94}\x94\x8c\x04rule\x94\x8c\x02cv\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c\x84/Users/bope/Documents/MutSig/scientafellow/packages/Projet/mutation_L1/python/m_results/manuscript/package/cucumber/snakemake/script\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/bope/Documents/MutSig/scientafellow/packages/Projet/mutation_L1/python/m_results/manuscript/package/cucumber/snakemake/script/ten_fold.py';
######## snakemake preamble end #########
#!/usr/bin/env python3
import pandas as pd
import numpy as np
from mutation_sim_train import main
from main_fixed_denovo import get_tri_context_fraction
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
