"""SEC 20220818 for assessing MPNN sequence outputs for designed coiled coils"""

"""Code from Deniz Akpinaroglu for aligning structures and getting RMSD """
import os
from pyrosetta import init
import pyrosetta.distributed.dask
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, OrResidueSelector
import pandas as pd

init()

####FUNCTIONS####

def calculate_rmsd(cc_pose, af_pose, file, h1start, h1end, h2start, h2end, min_name):
    all_rmsd_dict = {}
    helix1_all_rmsd_dict = {}
    helix1_rmsd_dict = {}
    helix2_all_rmsd_dict = {}
    helix2_rmsd_dict = {}

    #select helices only
    res_index_sel1 = ResidueIndexSelector()
    res_index_sel1.set_index_range(start=h1start, end=h1end)
    res_index_sel2 = ResidueIndexSelector()
    res_index_sel2.set_index_range(start=h2start, end=h2end)
    or_res_sel = OrResidueSelector()
    or_res_sel.add_residue_selector(res_index_sel1)
    or_res_sel.add_residue_selector(res_index_sel2)

    #select helix for ref only - add five to account for G4S loop
    res_index_sel2_ref = ResidueIndexSelector()
    res_index_sel2_ref.set_index_range(start=h2start+6, end=h2end+6)
    or_res_sel_ref = OrResidueSelector()
    or_res_sel_ref.add_residue_selector(res_index_sel1)
    or_res_sel_ref.add_residue_selector(res_index_sel2_ref)

    #bakcbone heavy atom RMSD for entire structure
    #align by both helices and get RMSD for whole structure
    all_rmsd_metric = RMSDMetric()
    all_rmsd_metric.set_comparison_pose(cc_pose)
    all_rmsd_metric.set_residue_selector_super(or_res_sel_ref)
    all_rmsd_metric.set_residue_selector_super_reference(or_res_sel)
    all_rmsd_metric.set_run_superimpose(True)
    all_rmsd_metric.set_residue_selector_reference(or_res_sel)
    all_rmsd_metric.set_residue_selector(or_res_sel_ref)
    all_rmsd_metric.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_atoms.
                                        rmsd_protein_bb_heavy_including_O)
    all_rmsd = all_rmsd_metric.calculate(af_pose)
    all_rmsd_dict[file] = all_rmsd

    #align by helix 1 and get RMSD for whole structure minus loop
    helix1_all_rmsd_metric = RMSDMetric()
    helix1_all_rmsd_metric.set_comparison_pose(cc_pose)
    helix1_all_rmsd_metric.set_residue_selector_super(res_index_sel1)
    helix1_all_rmsd_metric.set_run_superimpose(True)
    helix1_all_rmsd_metric.set_residue_selector_reference(or_res_sel)
    helix1_all_rmsd_metric.set_residue_selector(or_res_sel_ref)
    helix1_all_rmsd = helix1_all_rmsd_metric.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_atoms.
                                        rmsd_protein_bb_heavy_including_O)
    helix1_all_rmsd = helix1_all_rmsd_metric.calculate(af_pose)
    helix1_all_rmsd_dict[file] = helix1_all_rmsd

    #align by helix 1 and get RMSD for helix 1 only
    helix1_rmsd_metric = RMSDMetric()
    helix1_rmsd_metric.set_comparison_pose(cc_pose)
    helix1_rmsd_metric.set_residue_selector_super(res_index_sel1)
    helix1_rmsd_metric.set_run_superimpose(True)
    helix1_rmsd_metric.set_residue_selector(res_index_sel1)
    helix1_rmsd = helix1_rmsd_metric.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_atoms.
                                        rmsd_protein_bb_heavy_including_O)
    helix1_rmsd = helix1_rmsd_metric.calculate(af_pose)
    helix1_rmsd_dict[file] = helix1_rmsd

    #align by helix 2 and get RMSD for whole structure minus loop
    helix2_all_rmsd_metric = RMSDMetric()
    helix2_all_rmsd_metric.set_comparison_pose(cc_pose)
    helix2_all_rmsd_metric.set_residue_selector_super(res_index_sel2_ref)
    helix2_all_rmsd_metric.set_residue_selector_super_reference(res_index_sel2)
    helix2_all_rmsd_metric.set_run_superimpose(True)
    helix2_all_rmsd_metric.set_residue_selector_reference(or_res_sel)
    helix2_all_rmsd_metric.set_residue_selector(or_res_sel_ref)
    helix2_all_rmsd = helix2_all_rmsd_metric.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_atoms.
                                        rmsd_protein_bb_heavy_including_O)
    helix2_all_rmsd = helix2_all_rmsd_metric.calculate(af_pose)
    helix2_all_rmsd_dict[file] = helix2_all_rmsd

    #align by helix 2 and get RMSD for helix 2 only
    helix2_rmsd_metric = RMSDMetric()
    helix2_rmsd_metric.set_comparison_pose(cc_pose)
    helix2_rmsd_metric.set_residue_selector_super(res_index_sel2_ref)
    helix2_rmsd_metric.set_residue_selector_super_reference(res_index_sel2)
    helix2_rmsd_metric.set_run_superimpose(True)
    helix2_rmsd_metric.set_residue_selector_reference(res_index_sel2)
    helix2_rmsd_metric.set_residue_selector(res_index_sel2_ref)
    helix2_rmsd = helix2_rmsd_metric.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_atoms.
                                        rmsd_protein_bb_heavy_including_O)
    helix2_rmsd = helix2_rmsd_metric.calculate(af_pose)
    helix2_rmsd_dict[file] = helix2_rmsd


    all_rmsd_df = pd.DataFrame.from_dict(all_rmsd_dict, orient='index', columns=[f"{min_name}_all_rmsd"])
    all_rmsd_df.index.name = 'design'
    all_rmsd_df.reset_index(inplace=True)

    helix1_all_rmsd_df = pd.DataFrame.from_dict(helix1_all_rmsd_dict, orient='index', columns=[f"{min_name}_aligned_helix1_all_rmsd"])
    helix1_all_rmsd_df.index.name = 'design'
    helix1_all_rmsd_df.reset_index(inplace=True)

    helix1_rmsd_df = pd.DataFrame.from_dict(helix1_rmsd_dict, orient='index', columns=[f"{min_name}_aligned_helix1_rmsd"])
    helix1_rmsd_df.index.name = 'design'
    helix1_rmsd_df.reset_index(inplace=True)

    helix2_all_rmsd_df = pd.DataFrame.from_dict(helix2_all_rmsd_dict, orient='index', columns=[f"{min_name}_aligned_helix2_all_rmsd"])
    helix2_all_rmsd_df.index.name = 'design'
    helix2_all_rmsd_df.reset_index(inplace=True)

    helix2_rmsd_df = pd.DataFrame.from_dict(helix2_rmsd_dict, orient='index', columns=[f"{min_name}_aligned_helix2_rmsd"])
    helix2_rmsd_df.index.name = 'design'
    helix2_rmsd_df.reset_index(inplace=True)

    design_rmsd_df = pd.merge(all_rmsd_df, helix1_all_rmsd_df, on='design')
    design_rmsd_df = pd.merge(design_rmsd_df, helix1_rmsd_df, on='design')
    design_rmsd_df = pd.merge(design_rmsd_df, helix2_all_rmsd_df, on='design')
    design_rmsd_df = pd.merge(design_rmsd_df, helix2_rmsd_df, on='design')

    return design_rmsd_df

####
###USER DEFINED VARIABLES####
#get af output of all pdb files, align to
min1_bb_dir = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/min0_threads_dir'
colab_outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20230911_min0_52_ssd_bm01_loop_colabfold_array'
outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20231009_calc_af2_rmsd_with_loops_scripts'

helix_1_start = 1
helix_1_end = 28

helix_2_start = 29
helix_2_end = 70

####
colab_file_list = os.listdir(colab_outdir)
min1_file_list = os.listdir(min1_bb_dir)

# clfold pdb name: min1 then min 2: 20221010_30411_loop20221010_35251_loop_9_unrelaxed_rank_1_model_3.pdb
#input pdb name: 20221010_30411_loop.pdb

min1_dfs_to_concat = []

#220230321_16351_ALFA_52_9_unrelaxed_rank_005_alphafold2_multimer_v3_model_1_seed_000.pdb
#calculate RMSD to min1
for colab_file in colab_file_list:
    if colab_file.endswith('.pdb'):
        colab_pose = pyrosetta.pose_from_pdb(colab_outdir+'/'+colab_file) #af_pose
        print(len(colab_pose.sequence()))
        try:
            colab_pose = pyrosetta.pose_from_pdb(colab_outdir+'/'+colab_file) #af_pose
        except:
            pass
        colab_filename = colab_file.split('.')[0] #file
        min1_bb = colab_filename[0:22]
        min1_bb_filename = str(min1_bb+'.pdb')
        min1_pose = pyrosetta.pose_from_pdb(min1_bb_dir+'/'+min1_bb_filename) #cc_pose
        min1_rmsd_df = calculate_rmsd(cc_pose=min1_pose, af_pose=colab_pose, file=colab_filename,
                                h1start=helix_1_start, h1end=helix_1_end, h2start=helix_2_start, h2end=helix_2_end, min_name='min0')
        min1_dfs_to_concat.append(min1_rmsd_df)

print(min1_dfs_to_concat)
min1_df = pd.concat(min1_dfs_to_concat)

min1_df.to_csv(f"{outdir}/20231009_min0_52_df_all_outputs_bm01_loop.csv")

