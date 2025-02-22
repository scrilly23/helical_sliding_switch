import os
from frame2seq import Frame2seqRunner

####MAIN####

runner = Frame2seqRunner()

indir = os.getcwd()

#TODO: generalize, rn hardcoded for alfa tag at res 52/53
if '_52_' in indir:
    fixed_res_num = list(range(52, 65))
elif '_53_' in indir:
    fixed_res_num = list(range(53, 66))

for filename in os.listdir(indir):
    if filename.endswith('.pdb'):
        pdb_file = os.path.join(indir, filename)
        
        runner.design(pdb_file, chain_id='A', temperature=0.65, num_samples=10, omit_AA=['C', 'H'], fixed_positions=fixed_res_num, save_indiv_seqs=True)
