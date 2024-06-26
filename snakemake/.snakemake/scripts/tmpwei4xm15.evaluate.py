
######## snakemake preamble start (automatically inserted, do not edit) ########
#import sys; sys.path.extend(['/opt/homebrew/lib/python3.9/site-packages', '/Users/bope/Documents/MutSig/scientafellow/packages/Projet/mutation_L1/python/m_results/manuscript/mutational_starsign/snakemake/script']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xfb\x04\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c\x19results/data/1/test_m.csv\x94\x8c\x19results/data/1/11/0/s.csv\x94e}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x11\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x17)}\x94\x8c\x05_name\x94h\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x1dresults/data/1/11/0/e_new.csv\x94\x8c\x1eresults/data/1/11/0/logpmf.txt\x94e}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x04data\x94\x8c\x011\x94\x8c\x0211\x94\x8c\x010\x94e}\x94(h\r}\x94(\x8c\x07dataset\x94K\x00N\x86\x94\x8c\x04fold\x94K\x01N\x86\x94\x8c\x01k\x94K\x02N\x86\x94\x8c\x06Lambda\x94K\x03N\x86\x94uh\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94b\x8c\x07dataset\x94hE\x8c\x04fold\x94hFhOhG\x8c\x06Lambda\x94hHub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c0/var/folders/21/n1m80vz54z3djk8bv2d22z940000gn/T\x94e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bhgK\x01hiK\x01hkhdub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\x06config\x94}\x94\x8c\x04rule\x94\x8c\x08evaluate\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c|/Users/bope/Documents/MutSig/scientafellow/packages/Projet/mutation_L1/python/m_results/manuscript/mutational_starsign/snakemake/script\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/bope/Documents/MutSig/scientafellow/packages/Projet/mutation_L1/python/m_results/manuscript/mutational_starsign/snakemake/script/evaluate.py';
######## snakemake preamble end #########
#!/usr/bin/env python3
import pandas as pd
from starsigndna import refit
import numpy as np
snakemake = snakemake
test_m = pd.read_csv(snakemake.input[0], header=None, delimiter="\t")
test_m = test_m.iloc[:,1:]
test_m =  test_m.to_numpy().astype(float)
S = pd.read_csv(snakemake.input[1], header=None, delimiter=" ").to_numpy()
E,loss = refit(test_m, S , lambd=0)
#print(E)
np.savetxt(snakemake.output[0], E)
open(snakemake.output[1], "w").write(str(loss))
