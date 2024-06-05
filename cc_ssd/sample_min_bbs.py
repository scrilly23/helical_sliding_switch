
import pandas as pd
import numpy as np
import os
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
init()

from pyrosetta.rosetta.protocols.switches import GraftSwitchMover
import mdtraj as md

####FUNCTIONS####

def get_local_min_bbs(input_df, min_start, min_end):
    new_df = input_df.query('z1_offset >= @min_start and z1_offset <= @min_end', engine='python').copy()

    print('number of designs in local min: ' + str(list(new_df.shape)[0]))

    return new_df

def simple_thread_motif(input_pose, start_res, seq):
    thread_pose = pyrosetta.pose_from_pdb(input_pose)
    pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(f"""
    <SCOREFXNS>
    </SCOREFXNS>
    <MOVERS>
      <SimpleThreadingMover name="stm1"
                         pack_neighbors="false"
                         start_position= "{start_res}"
                         thread_sequence="{seq}"
                         pack_rounds="1"/>
    </MOVERS>
    """).get_mover("stm1").apply(thread_pose)
      
    return thread_pose

####

####USER DEFINED VARIABLES####
#evetually update these with argparse and os for making and managing dirs
#INPUTS
df = pd.read_csv('/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/20230622_0321_bbs_socket_cc_detected.csv')

bb_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20230321_bgs_sample_z1off_omega1_idealcc_valine'

fileheader = '20230627'

sf = get_fa_scorefxn()

min0_z_start = 8
min0_z_end = 11

min1_z_start = 13
min1_z_end = 16

min2_z_start = 18
min2_z_end = 21

tag_name = 'ALFA' #add in if statement so if not provided use motif sequence instead
seq2thread = 'SRLEEELRRRLTE'  #ALFA sequence
#important_residues = '3'  #ALFA important residues based on nanobody complex #SEC switching important residue to leucine only for burial
important_residues = '3,7'  #ALFA important residues for whole heptad burial
burial_cutoff = 6

#automate below by length of pose ###########
#length of whole cc - length of seq2thread - length of offset
pose_len = 70
seq2thread_len = len(seq2thread)
offset_len = min2_z_end - min1_z_start

#define limits of thread
start = pose_len - (seq2thread_len  + offset_len)
end = pose_len

#OUTPUTS
#need to make os if statement to make these dirs if they don't exist
outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins'
min0_thread_outdir = f"{outdir}/min0_threads_dir"
#min2_thread_outdir = f"{outdir}/min2_threads_dir"

#################################

#make dirs
outdirs = [min0_thread_outdir]

for dir in outdirs:
    if os.path.exists(dir) == False:
        os.makedirs(dir)

#for a cc with 7 heptad for a full turn, about 51 degrees from one heptad position to another (360/7)
#don't actually want helix to rotate that far
#for now, given the backbones I have up to 60 degree rotation, filtering out those up to 0.9 radians or 51.5 degrees
df = df.query('delta_omega1 < 0.9', engine='python').copy()

'''Sample backbones from user defined z1_offset ranges'''
#for now, not doing this by score, taking all backbones for which socket detects cc as defined by input df
min0_df = get_local_min_bbs(df, min0_z_start, min0_z_end)
min0_df['pdb_path'] = bb_indir + '/' + min0_df['structure']
min0_df.to_csv(f"{outdir}/{fileheader}_min0_bbs.csv", index=False)

'''Use GraftSwitchMover to place motif onto min0 backbones'''
#set up GraftSwitchMover
graftswitch = protocols.rosetta_scripts.XmlObjects.create_from_string(f"""
    <SCOREFXNS>
    </SCOREFXNS>
     <MOVERS>
        <GraftSwitchMover name="thread" sequence="{seq2thread}" important_residues="{important_residues}" start="{start}" burial_cutoff="{burial_cutoff}" end="{end}" pack_neighbors="true" pack_min="false"/>
    </MOVERS>
""").get_mover("thread")
graftswitch.score_function(sf)  # Set SF to REF15

min0_path_list = list(min0_df['pdb_path'])
min0_filename_list = list(min0_df['structure'])

for min0_bb_path, min0_file in zip(min0_path_list, min0_filename_list):
    min0_filename = min0_file.split('.pdb')[0]
    pose = pyrosetta.pose_from_pdb(min0_bb_path)
    graftswitch.apply(pose)

    #find first match solution
    first_match_position = pose.sequence().find(seq2thread) + 1
    if first_match_position == 0:
        continue
    elif first_match_position != 0:
        pose.dump_pdb(f"{min0_thread_outdir}/{min0_filename}_{tag_name}_{first_match_position}.pdb")

    # Repeat for any other solutions
    other_pose = graftswitch.get_additional_output()
    while(other_pose != None):
        match_position = other_pose.sequence().find(seq2thread) + 1
        if match_position == 0:
            continue
        elif match_position != 0:
            other_pose.dump_pdb(f"{min0_thread_outdir}/{min0_filename}_{tag_name}_{match_position}.pdb")
        other_pose = graftswitch.get_additional_output()




