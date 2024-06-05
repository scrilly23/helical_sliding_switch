
import pandas as pd
import os
import shutil

####USER DEFINED VARIABLES####
#evetually update these with argparse and os for making and managing dirs
#INPUTS
min2_df_no_clashes = pd.read_csv('/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_bbs_check_clashes/20230626_0321_min2_no_post_pack_bb_clashes.csv')

min0_header = 'min0'
min1_header = 'min1'
min2_header = 'min2'

min0_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/min0_threads_dir'
min1_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/min1_threads_dir'
min2_indir = '/wynton/home/kortemme/scrilly/helix_sliding/20230626_new_mins/min2_threads_dir'

#################################

src_dir = os.getcwd()

#get non redundant list of all thread positions with no clashes in min2 state
min2_df_no_clashes_temp = min2_df_no_clashes.copy()
min2_df_no_clashes_temp[['bb_id', 'thread_position']] = min2_df_no_clashes_temp['design'].str.split('_ALFA_', expand=True)
min2_df_no_clashes_temp[['thread_position', 'Nb_alfa_pdb']] = min2_df_no_clashes_temp['thread_position'].str.split('_', expand=True)
min2_df_no_clashes_temp = min2_df_no_clashes_temp.drop_duplicates(subset=['thread_position'])

min2_thread_solutions_list = min2_df_no_clashes_temp['thread_position'].tolist()
print(min2_thread_solutions_list)

#make min0, min1, min2 dirs  for each position
for i in min2_thread_solutions_list:
    dir_name = f'{min0_header}_{i}_threads'
    if os.path.exists(dir_name) == False:
        os.makedirs(dir_name)
    dir_name = f'{min1_header}_{i}_threads'
    if os.path.exists(dir_name) == False:
        os.makedirs(dir_name)
    dir_name = f'{min2_header}_{i}_threads'
    if os.path.exists(dir_name) == False:
        os.makedirs(dir_name)

#create min2 dir of threads without clashes for each thread position
design_id_list = min2_df_no_clashes.design.to_list()

new_design_id_list = []
for i in design_id_list:
    new_i = i[:-5]
    new_design_id_list.append(new_i)

pdb_string = '.pdb'
new_design_id_list = [x + pdb_string for x in new_design_id_list]

print(new_design_id_list)

#copy min2 bbs without clashes to new dir
min2_pdb_file_list = os.listdir(min2_indir)
for i in new_design_id_list:
    for j in min2_pdb_file_list:
        if i in j:
            for k in min2_thread_solutions_list:
                k_temp = f'_{k}.'
                if k_temp in j:
                    shutil.copy(f'{min2_indir}/{j}', f'{src_dir}/{min2_header}_{k}_threads')


#move files for min0 and min1 to new dirs by thread position
#need to clean this up and generalize
min0_pdb_file_list = os.listdir(min0_indir)
for j in min2_thread_solutions_list:
    j_temp = f'_{j}.'
    for file in min0_pdb_file_list:
        if j_temp in file:
            shutil.copy(f'{min0_indir}/{file}', f'{src_dir}/{min0_header}_{j}_threads')

min1_pdb_file_list = os.listdir(min1_indir)
for j in min2_thread_solutions_list:
    j_temp = f'_{j}.'
    for file in min1_pdb_file_list:
        if j_temp in file:
            shutil.copy(f'{min1_indir}/{file}', f'{src_dir}/{min1_header}_{j}_threads')






