import os
import pandas as pd

##### FUNCTIONS #####
#from Robert Alberstein MPNN_utils 2022
def MPNN_fasta_to_df(fname):
	"""
	Reads in an MPNN output fasta file and returns a pandas df with columns [T, sample, score, seq_recovery, sequence]
	"""
	all_data = []
	with open(fname, "r") as fin:
		# Store info for original sequence
		orig_info = fin.readline()[1:-1]  # omit starting '>' and newline character
		orig_score = orig_info.split()[1].split("=")[1][:-1]
		orig_seq = fin.readline()[:-1]  # omit newline character
		all_data.append([0.0, -1, float(orig_score), 1.0000, orig_seq])
		# Get info for rest of sequences
		for line in fin:
			if (line.startswith(">")):
				seq_info = line[1:-1].split(", ")
				all_data.append([
								float(seq_info[0].split("=")[-1]),
								int(seq_info[1].split("=")[-1]),
								float(seq_info[2].split("=")[-1]),
								float(seq_info[3].split("=")[-1])
								])
				all_data[-1].append(fin.readline()[:-1])

	df = pd.DataFrame(all_data, columns=["T", "sample", "score", "seq_recovery", "sequence"])
	df.index.name = 'Name'
	df.reset_index(inplace=True)

	return df


def write_df_to_fastas(df):
    """
	Write sequence DF to standard fasta file. Supply multi-row DFs to get combined FASTA files or single DF lines to do individual outputs.
	"""
    for index, row in df.iterrows():
        with open(f"{outdir}/{row['seq_id']}.fa", "w") as fout:
            one_chain = row['sequence'].split('/')[0]
            two_chain = row['sequence'].split('/')[1]
            fout.write(f">{row['seq_id']}\n{one_chain}SDPRKK{two_chain}\n")


####USER VARIABLES####
fasta_dir = '/wynton/home/kortemme/scrilly/helix_sliding/20230710_single_state_seqs/min0_52_mpnn_seqs/seqs'

outdir = '/wynton/home/kortemme/scrilly/helix_sliding/20231003_0321_min0_52_ssd_bm01_loop_sep_fa'

file_header = '20231003_0321_min0_52_ssd_bm01_loop'

####

os.makedirs(outdir, exist_ok=True)

#make and save df of all generated mpnn sequences

all_seqs_dfs_to_concat = []

for path, subdirs, files in os.walk(fasta_dir):
	print(path)
	for file in files:
		if file.endswith('.fa'):
			#print(file)
			filename = file.split('.')
			modelnames = filename[0]
			df = MPNN_fasta_to_df(path+'/'+file)
			df['model_names'] = modelnames
			#print(df)
			#print(df.shape)
			all_seqs_dfs_to_concat.append(df)
			df.to_csv(f"{fasta_dir}/{modelnames}_seqs.csv")


all_seqs_df = pd.concat(all_seqs_dfs_to_concat)
all_seqs_df['seq_id'] = all_seqs_df['model_names'] + '_' + all_seqs_df['Name'].astype(str)
all_seqs_df.to_csv(f"{outdir}/{file_header}_mpnn_ssd_all_seqs.csv")

#get rid of original seqs
new_seqs_df = all_seqs_df.query('seq_recovery < 1.0', engine='python').copy()
new_seqs_df.to_csv(f"{outdir}/{file_header}_mpnn_ssd_mpnn_seqs_only.csv")

#write to fasta for submision to Colabfold
write_df_to_fastas(new_seqs_df)
