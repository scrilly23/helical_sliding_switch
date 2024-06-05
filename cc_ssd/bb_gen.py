import os
import pandas as pd
from pyrosetta import *
pyrosetta.init()

####USER DEFINED VARIABLES####
sfxn = pyrosetta.rosetta.core.scoring.get_score_function() 

outdir = "/wynton/home/kortemme/scrilly/helix_sliding/20230321_bgs_sample_z1off_omega1_idealcc_valine"


##### FUNCTIONS #####
''' Generates an antiparallel helix pair of specified length from parametric equations '''
#set omega0 and r0 based on medians of naturally occurring cc's from Grigoryian and Degrado
#note that z1offset according to grigoryian and degrado occurs over a much smaller range in natural cc's than sampled here
def make_ap_helix_pair(num_res1, num_res2):
    pose = pyrosetta.Pose()
    pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(f"""
    <SCOREFXNS>
      <ScoreFunction name="sfxn1" weights="ref2015"/>
    </SCOREFXNS>
    <MOVERS>
      <BundleGridSampler name="bgs1"
                         max_samples="18000"
                         use_degrees= "1"
                         residue_name="VAL"
                         scorefxn="sfxn1"
                         r0="4.85"
                         omega0="-3.60"
                         dump_pdbs="1"
                         pdb_prefix= "20230321">
      <Helix
        crick_params_file="alpha_helix_100.crick_params"
        helix_length="{num_res1}"/>
      <Helix
        crick_params_file="alpha_helix_100.crick_params"
        helix_length="{num_res2}"
        delta_omega0="180"
        z1_offset_min="-5.00"
        z1_offset_max="25.00"
        z1_offset_samples="300"
        delta_omega1_min="-60.0"
        delta_omega1_max="60.0"
        delta_omega1_samples="60"
        invert="1"/>
    </BundleGridSampler>
    </MOVERS>
    """).get_mover("bgs1").apply(pose)
    return pose

'''For all pdb files in directory, pose, and score with ref2015
append pdb_file_name:score as key:value to dict'''
def pose_and_score_dir(directory):
    file_list = os.listdir(directory)
    score_dict = {}

    for file in file_list:
        if file.endswith('.pdb'):
            pose = pyrosetta.pose_from_pdb(directory+'/'+file)
            pdb_score = sfxn(pose)
            score_dict[file] = pdb_score

    return score_dict

    
    ########################################################
if (__name__ == "__main__"):
    # Make output dir if not current directory
    if (outdir != "."):
        os.makedirs(outdir, exist_ok=True)

    os.chdir(outdir)

    cc_pose = make_ap_helix_pair(28, 42)
    cc_pose.dump_pdb('original_cc_pose.pdb')

    os.chdir('..')

    cc_scores = pose_and_score_dir(outdir)

    file_list = os.listdir(outdir)
    z1off_dict = {}
    delta_omega1_dict = {}

    for file in file_list:
        if file.endswith('.pdb'):
            if file == 'original_cc_pose.pdb':
                continue
            else:
                #print(file)
                #print(type(file))
                with open(outdir+'/'+file, 'r') as f:

                    lines = f.readlines()

                    #here specifying lines with metrics sampled for helix 2 in a two helix cc
    				#change as needed depending on how BundleGridSampler was used
                    lines_helix2 = lines[42:52]

                    for line in lines_helix2:
                        if 'z1_offset' in line:
                            x = line.split(': ')
                            helix2_z1off = x[1]
                            z1off_dict[file] = helix2_z1off
                        if 'delta_omega1' in line:
                            y = line.split(': ')
                            helix2_delta_omega1 = y[1]
                            delta_omega1_dict[file] = helix2_delta_omega1

'''Create data frame of structures, z1offset, delta_omega1, and score'''
score_df = pd.DataFrame.from_dict(cc_scores, orient='index', columns=['score'])
score_df.index.name = 'structure'
score_df.reset_index(inplace=True)
#print(score_df)

z1off_df = pd.DataFrame.from_dict(z1off_dict, orient='index', columns=['z1_offset'])
z1off_df.index.name = 'structure'
z1off_df.reset_index(inplace=True)
#print(z1off_df)

delta_omega1_df = pd.DataFrame.from_dict(delta_omega1_dict, orient='index', columns=['delta_omega1'])
delta_omega1_df.index.name = 'structure'
delta_omega1_df.reset_index(inplace=True)
#print(delta_omega1_df)

score_z1off_df = pd.merge(score_df, z1off_df, on='structure')
#print(score_z1off_df)

score_z1off_delta_omega1_df = pd.merge(score_z1off_df, delta_omega1_df, on='structure')

#optional write to csv
score_z1off_delta_omega1_df.to_csv(outdir+'/20230321_bgs_sample_z1off_delta_omega1_valine.csv')