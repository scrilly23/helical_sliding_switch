import os
import pandas as pd
import argparse

##### FUNCTIONS #####
	
def fasta_to_df(fname, method='MPNN'):
    '''
    Converts a fasta file of sequences to a pandas dataframe.
    Adapted from Robert Alberstein's MPNN_utils 2022

    Parameters:
    fname (str): Path to fasta file. 
    
    method (str): Either 'MPNN' for ProteinMPNN or 'F2S' Frame2seq

    Returns:
    df: Returns a df with sequence and score info.
    '''
    #TODO: add for multistate
    all_data = []

    if method == 'MPNN':

        design_name = None
        with open(fname, "r") as fin:
            #get name of design backbone
            orig_info = fin.readline()[1:-1]
            design_name = orig_info.split()[0].split(", ")[0][1:-1]

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
        df['seq_id'] = design_name + '_' + df['Name'].astype(str)

    elif method == 'F2S':

        with open(fname, "r") as fin:
            design_name = fname.split('.fasta')[0].replace('A_seq', '') #TODO: fix hardcoded for chain A

            for line in fin:
                if (line.startswith(">")):
                    seq_info = line[1:-1].split(" ")
                    all_data.append([
					    seq_info[1].split("=")[-1],
					    seq_info[2].split("=")[-1], #TODO: seq_recovery-split on % and convert to float
					    float(seq_info[3].split("=")[-1]),
                        float(seq_info[4].split("=")[-1])
				    ])
                    all_data[-1].append(fin.readline()[:-1])

        df = pd.DataFrame(all_data, columns=["chain", "seq_recovery", "score", "T", "sequence"])
        df['seq_id'] = design_name

    return df


def write_df_to_fastas(df, loop_seq='', loop_cutpoint=None, single_chain_input=True):
    '''
    Writes out fasta files for each sequence in the dataframe.
    Adapted from Robert Alberstein's MPNN_utils 2022'
    
    Parameters:
    df: Dataframe containing sequence information.

    loop_seq (str): Optional. Sequence to add to the loop region of the designed sequence.

    loop_cutpoint (int): Optional. Index to cut the designed sequence at to add the loop sequence.

    single_chain_input (bool): Whether the input sequence is a single chain or not.

    Returns:
    new_seq_df: Dataframe with the new sequences with loops added if applicable.
    Directory in the current working directory with fasta files.
    '''
    cwd = os.getcwd()
    outdir = f'{cwd}/fastas_for_colabfold'
    os.mkdir(outdir)

    new_seq_dict = {}

    #TODO: make sure all combo options for design/predict are covered
    #need to add designed as multi chain, predicted as single chain no loop
    for index, row in df.iterrows():
        if single_chain_input:
            one_chain = row['sequence']
            #designed as single chain, predicted as single chain 
            #no sequence added
            if loop_cutpoint is None:
                print('oh no!')
                with open(f"{outdir}/{row['seq_id']}.fa", "w") as fout:
                    fout.write(f">{row['seq_id']}\n{one_chain}\n")
                new_seq_dict[row['seq_id']] = one_chain

            #designed as single chain, predicted as single chain 
            #sequence (eg. loop) added        
            else:
                with open(f"{outdir}/{row['seq_id']}.fa", "w") as fout:
                    fout.write(f">{row['seq_id']}\n{one_chain[:loop_cutpoint]}{loop_seq}{one_chain[loop_cutpoint:]}\n")
                new_seq_dict[row['seq_id']] = f"{one_chain[:loop_cutpoint]}{loop_seq}{one_chain[loop_cutpoint:]}"
        
        else:
            with open(f"{outdir}/{row['seq_id']}.fa", "w") as fout:
                one_chain = row['sequence'].split('/')[0]
                two_chain = row['sequence'].split('/')[1]
			    
                #designed as multichain, predicted as multichain
                if loop_seq == '':
                    fout.write(f">{row['seq_id']}\n{one_chain}:{two_chain}\n")
                    new_seq_dict[row['seq_id']] = f"{one_chain}:{two_chain}"
                
                #designed as multichain, predicted as single chain
                #loop sequence added
                else:
                    fout.write(f">{row['seq_id']}\n{one_chain}{loop_seq}{two_chain}\n")
                    new_seq_dict[row['seq_id']] = f"{one_chain}{loop_seq}{two_chain}"
    
    new_seq_df = pd.DataFrame.from_dict(new_seq_dict, orient='index', columns=['sequence'])
    new_seq_df.index.name = 'seq_id'
    new_seq_df.reset_index(inplace=True)

    return new_seq_df

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--method', type=str, default='MPNN')
    parser.add_argument('--loop', type=str, default='')
    parser.add_argument('--cutpoint', type=int, default=None)
    parser.add_argument('--single_chain', action='store_true')
    args = parser.parse_args()

    #process generated seqs
    indir = os.getcwd()
    files = os.listdir(indir)

    all_seqs_dfs_to_concat = []

    for file in files:
        if file.endswith(('.fa', '.fasta')):
            df = fasta_to_df(file, args.method)
            all_seqs_dfs_to_concat.append(df)
    
    all_seqs_df = pd.concat(all_seqs_dfs_to_concat)
    all_seqs_df.to_csv(f"{indir}/{args.method}_all_seqs_og.csv")

    new_seqs_df = write_df_to_fastas(all_seqs_df, loop_cutpoint=args.cutpoint, loop_seq=args.loop, single_chain_input=args.single_chain)
    new_seqs_df.to_csv(f"{indir}/{args.method}_all_seqs_for_colabfold.csv")

