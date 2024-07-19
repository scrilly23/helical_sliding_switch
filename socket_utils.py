import pandas as pd
import os


####FUNCTIONS####
def get_pdb_filename(input_fhread):
    
    pdb_identifier = None
    for line in input_fhread:
        if ('COILED COILS' in line) and ('result' in line):
            #print(line)
            pdb_identifier = line.split('.pdb')[0]

    return pdb_identifier

def df_from_dict(input_dict, column_name, index_name = 'design_id'):
    '''
    Converts a dictionary of saved values mapped to identifier to a df.
    
    Parameters:
    input_dict (dict): Dictionary with key corresponding to 
    design identifer and value some metric of interest.
    column_name (str): Name of data column for metric of interest.
     
    Returns:
    df: Returns df with one column for design_id and 
    the other for the metric of interest.
    
    '''
    df = pd.DataFrame.from_dict(input_dict, orient='index', columns=[column_name])
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': index_name}) 

    return df

def get_num_heptads(input_df):
    #TODO-should rework as list comprehension, rather than iterrows, will be faster
    '''
    Counts the number of heptads in each helix of a cc dimer from a socket output.

    Parameters:
    input_df (df): Data frame with columns containing columns with 
    helix 1 and helix 2 cc heptad register assignments from a socket output.

    Returns:
    df: Returns a df with a new column for each helix containing 
    the number of heptad repeats detected in the socket assignment for that helix. 
    '''

    h1_num_heptads_dict = {}
    h2_num_heptads_dict = {}

    for index, row in input_df.iterrows():
        h1_heptads = row['h1_reg'].count('abcdefg')
        h1_num_heptads_dict[row['design_id']] = h1_heptads
    
        h2_heptads = row['h2_reg'].count('abcdefg')
        h2_num_heptads_dict[row['design_id']] = h2_heptads

    h1_num_heptads_df = df_from_dict(h1_num_heptads_dict, column_name='h1_num_heptads')
    h2_num_heptads_df = df_from_dict(h2_num_heptads_dict, column_name='h2_num_heptads')

    new_df = input_df.merge(h1_num_heptads_df, how='left', on='design_id')
    new_df = new_df.merge(h2_num_heptads_df, how='left', on='design_id')

    return new_df

def get_num_ad(input_df):
    #TODO-should rework as list comprehension, rather than iterrows, will be faster
    '''
    Counts the number of a and d positions in each helix 
    of a cc dimer from a socket output.
    
    Added this in addition to counting heptads to account for 
    coiled coils with detected non-canonical heptad repeats.

    Parameters:
    input_df (df): Data frame with columns containing columns with 
    helix 1 and helix 2 cc heptad register assignments from a socket output.

    Returns:
    df: Returns a df with a new column for each helix containing 
    the number of a and d heptad positions detected 
    in the socket assignment for that helix. 
    '''

    h1_num_ad_dict = {}
    h2_num_ad_dict = {}

    for index, row in input_df.iterrows():
        h1_a = row['h1_reg'].count('a')
        h1_d = row['h1_reg'].count('d')
        h1_ad = h1_a + h1_d
        h1_num_ad_dict[row['design_id']] = h1_ad
    
        h2_a = row['h1_reg'].count('a')
        h2_d = row['h1_reg'].count('d')
        h2_ad = h2_a + h2_d
        h2_num_ad_dict[row['design_id']] = h2_ad

    h1_num_ad_df = df_from_dict(h1_num_ad_dict, column_name='h1_num_ad')
    h2_num_ad_df = df_from_dict(h2_num_ad_dict, column_name='h2_num_ad')

    new_df = input_df.merge(h1_num_ad_df, how='left', on='design_id')
    new_df = new_df.merge(h2_num_ad_df, how='left', on='design_id')

    return new_df

def heptad_to_upper(input_string, split_char = 'r '):
    split_1 = input_string.split(split_char)[0]
    split_2 = input_string.split(split_char)[1]
    split_2 = split_2.upper()
    
    new_string = split_1+split_char+split_2

    return new_string

def make_heptad_strings(reg_string, seq_string):
    a_list = []
    b_list = []
    c_list = []
    d_list = []
    e_list = []
    f_list = []
    g_list = []
    
    for index, char in enumerate(reg_string):
        if char == 'A':
            a_list.append(seq_string[index])
        elif char == 'B':
            b_list.append(seq_string[index])
        elif char == 'C':
            c_list.append(seq_string[index])
        elif char == 'D':
            d_list.append(seq_string[index])
        elif char == 'E':
            e_list.append(seq_string[index])
        elif char == 'F':
            f_list.append(seq_string[index])
        elif char == 'G':
            g_list.append(seq_string[index])
        else:
            continue
    
    a_string = ''.join(a_list)
    b_string = ''.join(b_list)
    c_string = ''.join(c_list)
    d_string = ''.join(d_list)
    e_string = ''.join(e_list)
    f_string = ''.join(f_list)
    g_string = ''.join(g_list)

    heptads_dict = {}
    heptads_dict['a_residues'] = a_string
    heptads_dict['b_residues'] = b_string
    heptads_dict['c_residues'] = c_string
    heptads_dict['d_residues'] = d_string
    heptads_dict['e_residues'] = e_string
    heptads_dict['f_residues'] = f_string
    heptads_dict['g_residues'] = g_string

    return heptads_dict


def get_heptad_res_identities(input_df):

    h1_dfs_to_concat = []
    h2_dfs_to_concat = []

    for index, row in input_df.iterrows():
        h1_reg = row['h1_reg']
        new_h1_reg = heptad_to_upper(h1_reg)

        h1_seq = row['h1_seq']

        h1_dict = make_heptad_strings(new_h1_reg, h1_seq)
        h1_df = pd.DataFrame([h1_dict])
        h1_df['h1_seq'] = h1_seq

        h1_dfs_to_concat.append(h1_df)

        h2_reg = row['h2_reg']
        new_h2_reg = heptad_to_upper(h2_reg)

        h2_seq = row['h2_seq']

        h2_dict = make_heptad_strings(new_h2_reg, h2_seq)
        h2_df = pd.DataFrame([h2_dict])
        h2_df['h1_seq'] = h1_seq

        h2_dfs_to_concat.append(h2_df)
    
    h1_dfs = pd.concat(h1_dfs_to_concat)
    h2_dfs = pd.concat(h2_dfs_to_concat)

    h1_h2_df = h1_dfs.merge(h2_dfs, on='h1_seq', suffixes=['_h1', '_h2'])

    new_df = input_df.merge(h1_h2_df, on='h1_seq')

    return new_df

def get_motif_res_in_heptad(input_df, motif = 'SRLEEELRRRLTE'):
    #should rework as list comprehension, rather than iterrows, will be faster
    '''
    Identifies the motif residues participating in heptad interaction of 
    helix2 of a cc dimer from a socket output.

    Parameters:
    input_df (df): Data frame with columns containing columns with 
    helix 1 and helix 2 cc heptad register assignments from a socket output.
    
    motif (str): A string of the binding motif of interest 
    threaded on to the coiled coil. 

    Returns:
    df: Returns a df with two new columns containing the motif residues 
    participating in the heptad interactions 
    and the number of motif residues participating in the heptad interactions. 
    '''
    motif_len = len(motif)

    seq_motif_res_dict = {}
    for index, row in input_df.iterrows():
        h2_seq = row['h2_seq']
        h2_reg = row['h2_reg']
        new_h2_reg = heptad_to_upper(h2_reg)

        motif_res_in_heptad = []

        motif_start = h2_seq.find(motif)

        if motif_start != -1:
            motif_end = motif_start + motif_len

            for i in range(motif_start, motif_end):
                if new_h2_reg[i] in ('A', 'B', 'C', 'D', 'E', 'F', 'G'):
                    motif_res_in_heptad.append(h2_seq[i]) 
                else:
                    continue

            motif_res_in_heptad_string = ''.join(motif_res_in_heptad)
            seq_motif_res_dict[h2_seq] = motif_res_in_heptad_string
        
        #handle case where full motif not detected in helical region
        elif motif_start == -1:
            sub_motif_start = len(h2_seq) - motif_len
            sub_motif_end = len(h2_seq)
            starting_sub_motif = motif[0:(sub_motif_end - sub_motif_start + 1)]

            j = 1
            for i in range(sub_motif_start, sub_motif_end):
                if h2_seq[i] != motif[0]:
                    pass
                elif h2_seq[i] == motif[0]: 
                    starting_sub_motif = motif[0:-j]
                    
                    if starting_sub_motif in h2_seq:
                        motif = starting_sub_motif
                    else:
                        pass
                
                j += 1
            
            motif_start = h2_seq.find(motif)
            motif_len = len(motif)

            motif_end = motif_start + motif_len

            for i in range(motif_start, motif_end):
                if new_h2_reg[i] in ('A', 'B', 'C', 'D', 'E', 'F', 'G'):
                    motif_res_in_heptad.append(h2_seq[i])
                else:
                    continue

                motif_res_in_heptad_string = ''.join(motif_res_in_heptad)
                seq_motif_res_dict[h2_seq] = motif_res_in_heptad_string

    seq_motif_res_df = df_from_dict(seq_motif_res_dict, index_name='h2_seq', column_name='motif_res_in_heptad')
    new_df = input_df.merge(seq_motif_res_df, how='left', on='h2_seq')
    new_df['num_motif_res_in_heptad'] = new_df['motif_res_in_heptad'].str.len()

    return new_df


if __name__ == '__main__':

#####USER DEFINED INPUTS#####
    indir='/wynton/home/kortemme/scrilly/helix_sliding/20240517_mpnn_msd_t_01_cf_array'
    #indir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240715_socket_utils_test'
    pdb_header = '13632_' #TODO-need to update this to be more general  
    outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20240517_mpnn_msd_t_01_cf_array'
    #outdir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20240715_socket_utils_test'
    file_header = 'MSD_13632_07144'
#####

    cc_dict = {}
    h1_seq_dict = {}
    h1_reg_dict = {}
    h2_seq_dict = {}
    h2_reg_dict = {}
    h1_non_canon_dict = {}
    h2_non_canon_dict = {}

    for path, subdirs, files in os.walk(indir):
        for file in files:
            if '.sh.o' in file:
                fh = open(path+'/'+file, encoding='utf-8', errors='ignore')
                fhread = fh.readlines()

                pdb_id = get_pdb_filename(fhread)

                for index, line in enumerate(fhread):
                    if (line.startswith(pdb_header)) and ('result' in line):
                        if 'NO COILED COILS' in line:
                            cc_dict[pdb_id] = 0
                        elif 'COILED COILS PRESENT' in line:
                            cc_dict[pdb_id] = 1
                        else:
                            cc_dict[pdb_id] = ''
                    if line.startswith('assigning heptad to helix 0'):
                        h1_seq = fhread[index+2]
                        h1_reg = fhread[index+3]
                        h1_seq_dict[pdb_id] = h1_seq
                        h1_reg_dict[pdb_id] = h1_reg

                        h1_non_canon = fhread[index+6]
                        h1_number_non_canon = h1_non_canon[10]
                        h1_non_canon_dict[pdb_id] = h1_number_non_canon
                    if line.startswith('assigning heptad to helix 1'):
                        h2_seq = fhread[index+2]
                        h2_reg = fhread[index+3]
                        h2_seq_dict[pdb_id] = h2_seq
                        h2_reg_dict[pdb_id] = h2_reg

                        h2_non_canon = fhread[index+6]
                        h2_number_non_canon = h2_non_canon[10]
                        h2_non_canon_dict[pdb_id] = h2_number_non_canon
                    else:
                        pass

                fh.close()


    list_of_dicts = [h1_seq_dict, h1_reg_dict, h2_seq_dict, h2_reg_dict, h1_non_canon_dict, h2_non_canon_dict]
    list_of_col_names = ['h1_seq', 'h1_reg', 'h2_seq', 'h2_reg', 'h1_non_canon_num_res', 'h2_non_canon_num_res']

    socket_call_df = df_from_dict(cc_dict, column_name='socket_call')

    for index, dict in enumerate(list_of_dicts):
        new_colname = list_of_col_names[index]
        temp_df = df_from_dict(dict, new_colname)

        socket_call_df = socket_call_df.merge(temp_df, how='left', on='design_id')

    socket_call_df.to_csv(f'{outdir}/{file_header}_socket_filtered.csv')

    og_num_seqs = socket_call_df.shape[0]

    #filter out designs without Socket2-detected kih
    socket_call_df = socket_call_df.query('socket_call == 1').copy()

    socket_pass_num_seqs = socket_call_df.shape[0]

    print(file_header)
    print(f'{socket_pass_num_seqs}/{og_num_seqs} seqs have kih interactions detected by Socket2.')

    df = get_num_heptads(socket_call_df)
    df = get_num_ad(socket_call_df)
    df = get_heptad_res_identities(socket_call_df)
    df = get_motif_res_in_heptad(socket_call_df)
    print(df)
    print(df.columns)


    df.to_csv(f'{outdir}/{file_header}_all_socket_outputs.csv')



