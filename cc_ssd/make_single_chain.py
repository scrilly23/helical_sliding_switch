from Bio.PDB import PDBParser, PDBIO
import os

#https://stackoverflow.com/questions/70246451/how-do-i-change-the-chain-name-of-a-pdb-file
#https://stackoverflow.com/questions/53094913/using-biopython-to-parse-a-pdb-file

#for given directories
#loop over pdbs and rename chain
#save into new dir

####USER DEFINED VARIABLES####
min0_52_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min0_52_threads'
min0_53_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min0_53_threads'

min2_52_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min2_52_threads'
min2_53_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20250220_min_bbs_for_f2s_ssd/min2_53_threads'

####

####FUNCTIONS####

def rename_chain(pdb_path, new_chain='A'):
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

####

####MAIN####

indir_list = [min0_52_indir, min0_53_indir, min2_52_indir, min2_53_indir]

for dir in indir_list:
    for filename in os.listdir(dir):
        if filename.endswith('.pdb'):
            single_chain_structure = rename_chain(os.path.join(dir, filename))
