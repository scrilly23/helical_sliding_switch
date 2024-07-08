import json
from pyrosetta.rosetta.protocols.moves import DsspMover
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, NotResidueSelector, OrResidueSelector
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core import *
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric

####FUNCTIONS####
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
    avg_plddt (flt)
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



#TODO-add def main, Karson rmsd, Dom rmsd






