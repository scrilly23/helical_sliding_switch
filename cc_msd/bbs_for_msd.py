
import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO
import mdtraj as md

####USER DEFINED VARIABLES####
min0_52_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min0_52_threads'
min0_53_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min0_53_threads'

min2_52_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min2_52_threads'
min2_53_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min2_53_threads'

outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20250415_bbs_for_msd'

designable_bbs_df_path = '/wynton/home/kortemme/scrilly/helix_sliding/20250415_bbs_for_msd/ssd_bbs_designable.csv'

####

####FUNCTIONS####

def rename_chain(pdb_path, new_chain='A'):
    '''
    Loop over pdbs in dir and rename chain
    
    Parameters:
    pdb_path (str): path to directory with pdb files.
    new_chain (str): chain id to rename to. Default is 'A'.
     
    Returns:
    Write out single chain pdbs to indir with same name.
    
    '''
    #TODO: generalize, add *args
    parser = PDBParser()
    io = PDBIO()

    structure = parser.get_structure('pdb_struct', pdb_path)

    for model in structure:
        for chain in model:
            chain.parent = None
            if chain == 'A':
                continue
            else:
                chain.id = new_chain

    io.set_structure(structure)

    io.save(pdb_path)


def merge_pdbs(path_to_pdb1, path_to_pdb2):
    '''
    Combine two pdbs into one pdb file for msd. 
    
    Parameters:
    path_to_pdb1 (str): path to first pdb file.
    path_to_pdb2 (str): path to second pdb file.
     
    Returns:
    Object with two pdbs combined. 
 
    '''
    p1 = md.load(str(path_to_pdb1))
    p2 = md.load(str(path_to_pdb2))
    top1, xyz1, time1, ul1, ua1 = p1.top, p1.xyz, p1.time, p1.unitcell_lengths, p1.unitcell_angles
    top2, xyz2, time2, ul2, ua2 = p2.top, p2.xyz, p2.time, p2.unitcell_lengths, p2.unitcell_angles
    top1_2 = top1.join(top2)

    translate_array = np.array([0,0,50])
    xyz3 = xyz2+translate_array
    xyz1_3 = np.concatenate((xyz1, xyz3), axis=1)
    p1_3 = md.Trajectory(xyz=xyz1_3, topology=top1_2)

    return p1_3


####

designable_bbs_df = pd.read_csv(designable_bbs_df_path)

indir_list = [min0_52_indir, min0_53_indir, min2_52_indir, min2_53_indir]

#get relevant pdbs and make single chain

for indir in indir_list:
    min_condition = indir[-15:-11]
    thread_position = indir[-10:-8]
    
    #make dir for each design set
    bb_outdir = os.path.join(outdir, f'{min_condition}_{thread_position}')
    
    if not os.path.exists(bb_outdir):
        os.makedirs(bb_outdir)

    thread_position = int(thread_position)

    if min_condition == 'min0':
        new_chain = 'A'
    elif min_condition == 'min2':
        new_chain = 'B'

    bbs = designable_bbs_df.query('min_condition == @min_condition & thread_position == @thread_position', engine='python')['bb_id']
    
    for filename in os.listdir(indir):
        if filename.endswith('.pdb'):
            for bb in bbs:
                if bb in filename:
                    #copy file to outdir
                    os.system(f'cp {indir}/{filename} {bb_outdir}/{filename}')

                    #rename chain
                    single_chain_structure = rename_chain(os.path.join(bb_outdir, filename), new_chain=new_chain)

#merge min0 and min2 pdbs to single file for msd

min0_52_pdb_files = os.listdir(f'{outdir}/min0_52')
min0_53_pdb_files = os.listdir(f'{outdir}/min0_53')

min2_52_pdb_files = os.listdir(f'{outdir}/min2_52')
min2_53_pdb_files = os.listdir(f'{outdir}/min2_53')

#make merged for thread 52
for file in min0_52_pdb_files:

    merged_bb_outdir = os.path.join(outdir, '52_threads_merged')
    
    if not os.path.exists(merged_bb_outdir):
        os.makedirs(merged_bb_outdir)

    if file.endswith('.pdb'):
        pdb1 = f'{outdir}/min0_52/{file}'
        pdb1_id = file[:-4]

        for file in min2_52_pdb_files:
            if file.endswith('.pdb'):
                pdb2 = f'{outdir}/min2_52/{file}'
                pdb2_id = file[:-4]

                merged_pdb = merge_pdbs(pdb1, pdb2)
                merged_pdb.save(f'{merged_bb_outdir}/{pdb1_id}_{pdb2_id}.pdb')


#make merged for thread 53
for file in min0_53_pdb_files:

    merged_bb_outdir = os.path.join(outdir, '53_threads_merged')
    
    if not os.path.exists(merged_bb_outdir):
        os.makedirs(merged_bb_outdir)

    if file.endswith('.pdb'):
        pdb1 = f'{outdir}/min0_53/{file}'
        pdb1_id = file[:-4]

        for file in min2_53_pdb_files:
            if file.endswith('.pdb'):
                pdb2 = f'{outdir}/min2_53/{file}'
                pdb2_id = file[:-4]

                merged_pdb = merge_pdbs(pdb1, pdb2)
                merged_pdb.save(f'{merged_bb_outdir}/{pdb1_id}_{pdb2_id}.pdb')

