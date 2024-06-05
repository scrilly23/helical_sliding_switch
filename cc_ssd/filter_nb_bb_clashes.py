import os
from pyrosetta import *
from pyrosetta import init
import pyrosetta.distributed.dask
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.simple_metrics.metrics import ResidueSummaryMetric, summary_type
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import PerResidueClashMetric
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, TrueResidueSelector
from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
from pyrosetta.rosetta.protocols import minimization_packing as pack_min

import pandas as pd

init('-ex1 -ex2')

#append 612g pose to each pose in min2 alfa thread dir-done

#superimpose_pose-done

#clash check on nanobody

#filter out min2 states with clashes

#get bb only clashes-done
#delete alfa tag portion from nb alfa structure after aligning but before calculating clashes-done
#switch to report clashes from cc rather than nb-done

####FUNCTIONS####

def superimpose_poses(pose1, pose1_start, pose1_end, pose2_start, pose2_end, pose2):
    pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(f"""
    <MOVERS>
      <Superimpose name="super_pose"
                         ref_start="{pose1_start}"
                         ref_end="{pose1_end}"
                         ref_pose="{pose1}"
                         target_start="{pose2_start}"
                         target_end="{pose2_end}"
                         CA_only="0"/>
    </MOVERS>
    """).get_mover("super_pose").apply(pose2)
    return pose2

#below adapted from https://github.com/proleu/hinge_paper/blob/main/crispy_shifty/protocols/states.py
def bb_clash_check(pose: Pose):
    all_gly = pose.clone()
    true_sel = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    true_x = true_sel.apply(all_gly)
    pyrosetta.rosetta.protocols.toolbox.pose_manipulation.repack_these_residues(true_x, all_gly, sfxn, False, 'G')
    clash = clash_metric.calculate(all_gly)
    clash_res = str(clash) 
    clash_sum = summary_metric.calculate(all_gly)

    return clash_res, clash_sum


###USER DEFINED VARIABLES####
#get af output of all pdb files, align to
nb_pdb_path = '/wynton/home/kortemme/scrilly/helix_sliding/6i2g_clean_one_chain.pdb'
min2_bb_dir = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/min2_threads_dir'
no_clash_pdb_outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_bbs_check_clashes'

sfxn = pyrosetta.get_fa_scorefxn()
####
min2_file_list = os.listdir(min2_bb_dir)

# clfold pdb name: min1 then min 2: 20221010_30411_loop20221010_35251_loop_9_unrelaxed_rank_1_model_3.pdb
#input pdb name: 20221010_30411_loop.pdb

pose_6i2g = pyrosetta.pose_from_pdb(nb_pdb_path) 
seq_6i2g = pose_6i2g.sequence()
pose_6i2g_alfa_start = seq_6i2g.find('SRLEEELRRRLTE') + 1
pose_6i2g_alfa_end = pose_6i2g_alfa_start + 12
print(pose_6i2g_alfa_start)
print(pose_6i2g_alfa_end)
print(len(seq_6i2g))
print(pose_6i2g.size())

#align 6i2g to cc by alfa tag residues
for min2_file in min2_file_list:
    min2_path = f"{min2_bb_dir}/{min2_file}"
    if min2_file.endswith('.pdb'):
        min2_pose = pyrosetta.pose_from_pdb(f"{min2_bb_dir}/{min2_file}")
        min2_filename = min2_file.split('.pdb')[0] #file
        
        min2_pose_seq = min2_pose.sequence()
        min2_pose_cc_end = len(min2_pose_seq)
        min2_pose_alfa_start = min2_pose_seq.find('SRLEEELRRRLTE') + 1
        min2_pose_alfa_end = min2_pose_alfa_start + 12
        
        new_pose = superimpose_poses(min2_path, min2_pose_alfa_start, min2_pose_alfa_end, pose_6i2g_alfa_start, pose_6i2g_alfa_end, pose_6i2g) #for some reason could not open 6i2g if provided as ref pose
        #get rid of alfa tag from 6i2g structure
        delete_mover = pyrosetta.rosetta.protocols.grafting.simple_movers.DeleteRegionMover()
        delete_mover.region(f"{pose_6i2g_alfa_start}", f"{pose_6i2g_alfa_end}")
        new_pose_no_alfa = new_pose.clone()
        delete_mover.apply(new_pose_no_alfa)
        pyrosetta.rosetta.core.pose.append_pose_to_pose(min2_pose, new_pose_no_alfa)

        min2_pose.dump_pdb(f"{no_clash_pdb_outdir}/{min2_filename}_6i2g.pdb")

combo_pdb_file_list = os.listdir(no_clash_pdb_outdir)

pre_pack_clash_dict = {}
pre_pack_clash_residues_dict = {}

pre_pack_bb_clash_dict = {}
pre_pack_bb_clash_residues_dict = {}

post_pack_clash_dict = {}
post_pack_clash_residues_dict = {}

post_pack_bb_clash_dict = {}
post_pack_bb_clash_residues_dict = {}

for combo_pdb_file in combo_pdb_file_list:
    if combo_pdb_file.endswith('6i2g.pdb'):
        combo_filename_split = combo_pdb_file.split('.pdb')
        combo_filename = combo_filename_split[0]
        combo_pose = pyrosetta.pose_from_pdb(f"{no_clash_pdb_outdir}/{combo_pdb_file}")
        combo_pose_seq = combo_pose.sequence()
        #search for and select nb alfa residues in cc multimer prediction
        cc_pose_nb_start = combo_pose_seq.find('VQLQESGGGLVQPGGSLRLSCTASGVTISALNAMAMGWYRQAPGERRVMVAAVSERGNAMYRESVQGRFTVTRDFTNKMVSLQMDNLKPEDTAVYYCHVLEDRVDSFHDYWGQGTQVTVSS') + 1
        cc_pose_nb_end = cc_pose_nb_start + 119
        print(cc_pose_nb_start)
        print(cc_pose_nb_end)
        nb_alfa_in_cc_sel = ResidueIndexSelector()
        nb_alfa_in_cc_sel.set_index_range(start=cc_pose_nb_start, end=cc_pose_nb_end)

        cc_sel = ResidueIndexSelector()
        cc_sel.set_index_range(start=1, end=min2_pose_cc_end)

        clash_metric = PerResidueClashMetric()
        clash_metric.set_residue_selector(cc_sel)
        clash_metric.set_secondary_residue_selector(nb_alfa_in_cc_sel)
        pre_pack_clash = clash_metric.calculate(combo_pose)

        pre_pack_clash_res = str(pre_pack_clash) #fine for now, but figure out how to unpack in a better way later
        pre_pack_clash_residues_dict[combo_filename] = pre_pack_clash_res

        summary_metric = ResidueSummaryMetric()
        summary_metric.set_metric(clash_metric)
        summary_metric.set_action(summary_type.sum)
        pre_pack_clash_sum = summary_metric.calculate(combo_pose)
        pre_pack_clash_dict[combo_filename] = pre_pack_clash_sum

        pre_pack_bb_clash_res, pre_pack_bb_clash_sum = bb_clash_check(combo_pose)
        pre_pack_bb_clash_residues_dict[combo_filename] = pre_pack_bb_clash_res
        pre_pack_bb_clash_dict[combo_filename] = pre_pack_bb_clash_sum

        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(operation.RestrictToRepacking())
        packer = pack_min.PackRotamersMover()
        packer.task_factory(tf)
        packer.apply(combo_pose)

        post_pack_clash = clash_metric.calculate(combo_pose)
        
        post_pack_clash_res = str(post_pack_clash)
        post_pack_clash_residues_dict[combo_filename] = post_pack_clash_res

        post_pack_clash_sum = summary_metric.calculate(combo_pose)
        post_pack_clash_dict[combo_filename] = post_pack_clash_sum

        post_pack_bb_clash_res, post_pack_bb_clash_sum = bb_clash_check(combo_pose)
        post_pack_bb_clash_residues_dict[combo_filename] = post_pack_bb_clash_res
        post_pack_bb_clash_dict[combo_filename] = post_pack_bb_clash_sum

#print(pre_pack_clash_dict)
#print(post_pack_clash_dict)
pre_pack_df = pd.DataFrame.from_dict(pre_pack_clash_dict, orient='index', columns=['pre_pack_clashes'])
pre_pack_df.index.name = 'design'
pre_pack_df.reset_index(inplace=True)

pre_pack_res_df = pd.DataFrame.from_dict(pre_pack_clash_residues_dict, orient='index', columns=['pre_pack_clashes_residues'])
pre_pack_res_df.index.name = 'design'
pre_pack_res_df.reset_index(inplace=True)

post_pack_df = pd.DataFrame.from_dict(post_pack_clash_dict, orient='index', columns=['post_pack_clashes'])
post_pack_df.index.name = 'design'
post_pack_df.reset_index(inplace=True)

post_pack_res_df = pd.DataFrame.from_dict(post_pack_clash_residues_dict, orient='index', columns=['post_pack_clashes_residues'])
post_pack_res_df.index.name = 'design'
post_pack_res_df.reset_index(inplace=True)

pre_pack_bb_df = pd.DataFrame.from_dict(pre_pack_bb_clash_dict, orient='index', columns=['pre_pack_bb_clashes'])
pre_pack_bb_df.index.name = 'design'
pre_pack_bb_df.reset_index(inplace=True)

pre_pack_bb_res_df = pd.DataFrame.from_dict(pre_pack_bb_clash_residues_dict, orient='index', columns=['pre_pack_bb_clashes_residues'])
pre_pack_bb_res_df.index.name = 'design'
pre_pack_bb_res_df.reset_index(inplace=True)

post_pack_bb_df = pd.DataFrame.from_dict(post_pack_bb_clash_dict, orient='index', columns=['post_pack_bb_clashes'])
post_pack_bb_df.index.name = 'design'
post_pack_bb_df.reset_index(inplace=True)

post_pack_bb_res_df = pd.DataFrame.from_dict(post_pack_bb_clash_residues_dict, orient='index', columns=['post_pack_bb_clashes_residues'])
post_pack_bb_res_df.index.name = 'design'
post_pack_bb_res_df.reset_index(inplace=True)

all_clashes_df = pre_pack_df.merge(post_pack_df, on='design', how='inner')
all_clashes_df = all_clashes_df.merge(pre_pack_res_df, on='design', how='inner')
all_clashes_df = all_clashes_df.merge(post_pack_res_df, on='design', how='inner')

all_clashes_df = all_clashes_df.merge(pre_pack_bb_df, on='design', how='inner')
all_clashes_df = all_clashes_df.merge(post_pack_bb_df, on='design', how='inner')
all_clashes_df = all_clashes_df.merge(pre_pack_bb_res_df, on='design', how='inner')
all_clashes_df = all_clashes_df.merge(post_pack_bb_res_df, on='design', how='inner')
#print(all_clashes_df['pre_pack_clashes'].min())
#print(all_clashes_df['post_pack_clashes'].min())
all_clashes_df.to_csv(f"{no_clash_pdb_outdir}/20230626_0321_min2_nb_clashes_ex_rot_df.csv")


