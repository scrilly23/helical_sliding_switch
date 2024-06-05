import pandas as pd

####USER DEFINED VARIABLES####
input_file_path = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20230621_socket_test/submit_20230621_socket_0326_bbs_v3_resaved.txt'
pdb_header = '20230321_'
outdir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20230621_socket_test'
file_header = '20230622'
####

fh = open(input_file_path)

fhread = fh.readlines()

cc_dict = {}

for line in fhread:
    if line.startswith(pdb_header):
        pdb_id = line.split('.pdb')[0]
        if 'NO COILED COILS' in line:
            cc_dict[pdb_id] = 0
        elif 'COILED COILS PRESENT' in line:
            cc_dict[pdb_id] = 1
        else:
            cc_dict[pdb_id] = ''

fh.close()

df = pd.DataFrame.from_dict(cc_dict, orient='index', columns=['socket_call'])
df.reset_index(inplace=True)
df = df.rename(columns = {'index':'design_id'})
print(df.shape)
df.to_csv(f'{outdir}/{file_header}_0321_bbs_socket_filtered.csv')

    

