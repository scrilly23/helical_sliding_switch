import json
import os
import pandas as pd
from pyrosetta import init
import pyrosetta.distributed.dask
from pyrosetta.rosetta.protocols.moves import DsspMover
from pyrosetta.rosetta.protocols.fold_from_loops.movers import AlignByResidueSelectorMover
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, NotResidueSelector, OrResidueSelector
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core import *
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric

pyrosetta.init()

####FUNCTIONS####
def df_from_dict(input_dict, column_name, index_name = 'design_id'):
    '''
    Converts a dictionary of saved values mapped to identifier to a df.
    
    Parameters:
    input_dict (dict): Dictionary with key corresponding to design identifer and value some metric of interest.
    column_name (str): Name of data column for metric of interest.
     
    Returns:
    df: Returns df with one column for design_id and the other for the metric of interest.
    
    '''
    df = pd.DataFrame.from_dict(input_dict, orient='index', columns=[column_name])
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': index_name}) 

    return df

def get_loop_resn(pose):
    #TODO-handle exception if there are no loops
    '''get loop residue numbers of pose and return as list

    Parameters:
    pose (pose): pose to get loop residue numbers from

    Returns:
    loop_resn (lst of ints): list of residue numbers of loops
    loop_resn_rosetta (lst of ints): list of residue numbers of loops in pose in pose numbering
    '''
    dssp_mover = DsspMover()
    dssp_mover.apply(pose)
    ss = pose.secstruct()

    #get index of resis that are in loops
    loop_resn = [i for i, x in enumerate(ss) if x == 'L']
    loop_resn_rosetta = [i+1 for i, x in enumerate(ss) if x == 'L']

    return loop_resn, loop_resn_rosetta


def calculate_rmsd_no_loops(rosetta_pose, af_pose):
    '''aligns two poses and calculates RMSD over non-loop residues

    Parameters:
    rosetta_pose (pose): rosetta model pose
    af_pose (pose): af2 predicted structure

    Returns:
    all_rmsd (flt): bb heavy atom rmsd over non-loop residues
    '''

    #get loop residue numbers from rosetta pose
    loop_resn, loop_resn_rosetta = get_loop_resn(rosetta_pose)

    #select non-loop residues of rosetta pose
    loop_selector = ResidueIndexSelector()
    for i in loop_resn_rosetta:
        loop_selector.append_index(i)

    non_loop_selector = NotResidueSelector(loop_selector)

    #align and rmsd over non-loop residues
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_comparison_pose(rosetta_pose)
    rmsd_metric.set_residue_selector_super(non_loop_selector)
    rmsd_metric.set_residue_selector_super_reference(non_loop_selector)
    rmsd_metric.set_run_superimpose(True)
    rmsd_metric.set_residue_selector_reference(non_loop_selector)
    rmsd_metric.set_residue_selector(non_loop_selector)
    rmsd_metric.set_rmsd_type(rmsd_atoms.rmsd_protein_bb_heavy_including_O)
    all_rmsd = rmsd_metric.calculate(af_pose)
    
    return all_rmsd


def calculate_cc_rmsd(rosetta_pose, af_pose):
    '''aligns two cc poses and calculates RMSD over non-loop residues

    Parameters:
    rosetta_pose (pose): rosetta model pose
    af_pose (pose): af2 predicted structure

    Returns:
    all_rmsd (flt): bb heavy atom rmsd over non-loop residues
    '''
    #calculate loop length by difference in length between poses
    loop_len = len(af_pose.sequence()) - len(rosetta_pose.sequence())
    
    #get resn rosetta pose helices stop and start
        #TODO-generalize this for any number of helices
    jump_resn = pose.chain_end_res(rosetta_pose)
    h1_start = 1
    h1_end = jump_resn[1]
    h2_start = jump_resn[1]+1
    h2_end = jump_resn[2]

    #select helix residues from rosetta cc pose
    cc_ref_res_sel = ResidueIndexSelector()
    cc_ref_res_sel.set_index_range(start=h1_start, end=h2_end)
    
    #select helix residues for af_pose w/out loop
    af_sel1 = ResidueIndexSelector()
    af_sel1.set_index_range(start=h1_start, end=h1_end)
    af_sel2 = ResidueIndexSelector()
    af_sel2.set_index_range(start=h2_start+loop_len, end=h2_end+loop_len)
    or_res_sel = OrResidueSelector()
    or_res_sel.add_residue_selector(af_sel1)
    or_res_sel.add_residue_selector(af_sel2)

    #align and rmsd over non-loop residues
    all_rmsd_metric = RMSDMetric()
    all_rmsd_metric.set_comparison_pose(rosetta_pose)
    all_rmsd_metric.set_residue_selector_super(or_res_sel)
    all_rmsd_metric.set_residue_selector_super_reference(cc_ref_res_sel)
    all_rmsd_metric.set_run_superimpose(True)
    all_rmsd_metric.set_residue_selector_reference(cc_ref_res_sel)
    all_rmsd_metric.set_residue_selector(or_res_sel)
    all_rmsd_metric.set_rmsd_type(rmsd_atoms.rmsd_protein_bb_heavy_including_O)
    all_rmsd = all_rmsd_metric.calculate(af_pose)

    return all_rmsd


def calculate_average(lst): 
    return sum(lst) / len(lst) 


def calculate_avg_plddt(json_file_path, loop_resn=None):
    '''get average plddt from colabfold output json.
    
    Parameters:
    json_file_path (str): filepath to the json file containing plddt values per residue. 
    loop_resn (lst of ints): residue numbers zero indexed of residues to ignore for calculating avg plddt, eg loops. 

    Returns:
    avg_plddt (flt)
    '''

    f = open(json_file_path,)
    data = json.load(f)
    plddt_vals = data['plddt']

    #plddt without excluding any regions
    if loop_resn == None:
        avg_plddt = calculate_average(plddt_vals)
    else:
        #plddt excluding loops, and/or other regions
        for i in sorted(loop_resn, reverse=True):
            del plddt_vals[i]
        avg_plddt = calculate_average(plddt_vals)
    
    return avg_plddt


def calculate_avg_pae(json_file_path, loop_resn=None):
    '''get average plddt from colabfold output json.
    
    Parameters:
    json_file_path (str): filepath to the json file containing pae values per residue. 
    loop_resn (lst of ints): residue numbers zero indexed of residues to ignore for calculating avg plddt, eg loops. 

    Returns:
    avg_pae (flt)
    '''

    f = open(json_file_path,)
    data = json.load(f)
    pae_vals = data['pae']

    #pae without excluding any regions
    if loop_resn == None:
        avg_pae_vals = []
        for paes in pae_vals:
            avg_pae = calculate_average(paes)
            avg_pae_vals.append(avg_pae)
        total_avg_pae = calculate_average(avg_pae_vals)
    else:
        #pae excluding loops, and/or other regions
        for i in sorted(loop_resn, reverse=True):
            del pae_vals[i]
        for paes in pae_vals:
            avg_pae = calculate_average(paes)
            avg_pae_vals.append(avg_pae)
        total_avg_pae = calculate_average(avg_pae_vals)
    
    return total_avg_pae

def align_by_residues(input_pose, nmr_pose): #TODO-generalize and make a cc version
    """
    @dwgrisingher
    Aligns two poses by residue selector over the whole protein
    reason I am doing it this way incase I have to align by specefic residues in the
    future I can just add them here
    """
     #selecting residues in LUCS pose
    original_res_1 = ResidueIndexSelector("1-28")
    original_res_2 = ResidueIndexSelector("29-70")

    original_res = OrResidueSelector()
    original_res.add_residue_selector(original_res_1)
    original_res.add_residue_selector(original_res_2)

    #selecting residues in NMR pose
    AF_res_1 = ResidueIndexSelector("1-28")
    AF_res_2 = ResidueIndexSelector("34-75")

    AF_res = OrResidueSelector()
    AF_res.add_residue_selector(AF_res_1)
    AF_res.add_residue_selector(AF_res_2)

    #print(f"Aligning mobile pose residues: {list(select.get_residue_set_from_subset(mobile_selector.apply(mobile_pose)))}")
    #print(f"to reference pose residues: {list(select.get_residue_set_from_subset(ref_selector.apply(ref_pose)))}")
    abr = AlignByResidueSelectorMover()
    abr.query_selector(AF_res)
    abr.reference_pose(input_pose)
    abr.reference_selector(original_res)
    abr.apply(nmr_pose)

def calculate_RMSDs_on_beta_only_helix(input_pose, AF_pose): #TODO-generalize and make a cc version
    """
    @dwgrisingher
    Takes two already-aligned poses and returns bb RMSDs
    """

    #selecting residues in LUCS pose
    original_res_1 = ResidueIndexSelector("1-28")
    original_res_2 = ResidueIndexSelector("29-70")

    original_res = OrResidueSelector()
    original_res.add_residue_selector(original_res_1)
    original_res.add_residue_selector(original_res_2)

    #selecting residues in NMR pose
    AF_res_1 = ResidueIndexSelector("1-28")
    AF_res_2 = ResidueIndexSelector("34-75")

    AF_res = OrResidueSelector()
    AF_res.add_residue_selector(AF_res_1)
    AF_res.add_residue_selector(AF_res_2)

    # Set up the RMSDMetric
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_run_superimpose(False)
    rmsd_metric.set_comparison_pose(input_pose)
    rmsd_metric.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_atoms.rmsd_protein_bb_heavy_including_O)

    # Calculate RMSD of whole structure
    rmsd_metric.set_residue_selector(AF_res)
    rmsd_metric.set_residue_selector_reference(original_res)
    rmsd = rmsd_metric.calculate(AF_pose)

    return rmsd

#TODO-add def main, Karson rmsd

if __name__ == "__main__":

    ###USER DEFINED VARIABLES####
    #get af output of all pdb files, align to
    min0_path = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/min0_threads_dir/20230321_13632_ALFA_52.pdb'
    min2_path = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/min2_threads_dir/20230321_07144_ALFA_52.pdb'
    colab_outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20240517_mpnn_msd_t_01_cf_array'
    outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20240517_mpnn_msd_t_01_cf_array'
    ####

    #load starting state poses
    min0_pose = pyrosetta.pose_from_pdb(min0_path)
    min2_pose = pyrosetta.pose_from_pdb(min2_path)

    colab_file_list = os.listdir(colab_outdir)

    min0_rmsd_dict = {}
    min2_rmsd_dict = {}
    avg_plddt_dict = {}
    avg_pae_dict = {}

    #calculate RMSD to min0 and min2 states
    for colab_file in colab_file_list:

        if colab_file.endswith('.pdb'):
            #get rmsd
            colab_pose = pyrosetta.pose_from_pdb(colab_outdir+'/'+colab_file) #af_pose

            colab_filename = colab_file.split('.')[0] #file

            rmsd_to_min0 = calculate_cc_rmsd(min0_pose, colab_pose)
            rmsd_to_min2 = calculate_cc_rmsd(min2_pose, colab_pose)
            
            min0_rmsd_dict[colab_filename] = rmsd_to_min0
            min2_rmsd_dict[colab_filename] = rmsd_to_min2

            #get avg plddt across non-loop residues
            loop_residues, rosetta_loop_residues = get_loop_resn(colab_pose)

            json_file = colab_file.replace('.pdb', '.json')
            json_file = json_file.replace('unrelaxed', 'scores')

            avg_plddt = calculate_avg_plddt(colab_outdir+'/'+json_file, loop_residues)

            avg_plddt_dict[colab_filename] = avg_plddt

            #get avg pae across non-loop residues
            pae_file = colab_file.split('unrelaxed')[0]
            pae_file = pae_file + 'predicted_aligned_error_v1.json'

            avg_pae = calculate_avg_pae(colab_outdir+'/'+pae_file, loop_residues)

            avg_pae_dict[colab_filename] = avg_pae
    
    min0_rmsd_df = df_from_dict(input_dict=min0_rmsd_dict, column_name='min0_all_rmsd_no_loop')
    min2_rmsd_df = df_from_dict(input_dict=min2_rmsd_dict, column_name='min2_all_rmsd_no_loop')
    avg_plddt_df = df_from_dict(input_dict=avg_plddt_dict, column_name='avg_plddt_no_loop')
    avg_pae_df = df_from_dict(input_dict=avg_pae_dict, column_name='avg_pae_no_loop')

    all_df = pd.merge(min0_rmsd_df, min2_rmsd_df, on='design_id')
    all_df = pd.merge(all_df, avg_plddt_df, on='design_id')
    all_df = pd.merge(all_df, avg_pae_df, on='design_id')

    all_df.to_csv(f'{outdir}/af2_metrics.csv')

#pdb file: 13632_ALFA_52_07144_ALFA_52_t_01_9_unrelaxed_rank_005_alphafold2_multimer_v3_model_5_seed_000.pdb
#pae file: 13632_ALFA_52_07144_ALFA_52_t_01_999_predicted_aligned_error_v1.json
#plddt file: 13632_ALFA_52_07144_ALFA_52_t_01_99_scores_rank_001_alphafold2_multimer_v3_model_3_seed_000.json





