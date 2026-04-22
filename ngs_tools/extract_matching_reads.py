import gzip
import pandas as pd
import time
import os
import glob
import sys
from argparse import ArgumentParser
from typing import Iterator, Tuple

def parse_fastq_gz(filepath: str) -> Iterator[Tuple[str, str, str, str]]:
    """
    Parse a FASTQ.gz file and yield records as tuples.
    
    Args:
        filepath: Path to the FASTQ.gz file
        
    Yields:
        Tuple of (header, sequence, plus_line, quality) for each record
    @Claude
    """
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        while True:
            # Read 4 lines at a time (FASTQ format)
            header = f.readline().strip()
            if not header:  # End of file
                break
                
            sequence = f.readline().strip()
            plus_line = f.readline().strip()
            quality = f.readline().strip()
            
            # Validate FASTQ format
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ header: {header}")
            if not plus_line.startswith('+'):
                raise ValueError(f"Invalid FASTQ plus line: {plus_line}")
                
            yield header, sequence, plus_line, quality

def fastq_gz_to_dataframe(filepath: str, max_records: int = None) -> pd.DataFrame:
    """
    Convert a FASTQ.gz file to a pandas DataFrame.
    
    Args:
        filepath: Path to the FASTQ.gz file
        max_records: Optional limit on number of records to read
        
    Returns:
        pandas DataFrame with columns: sequence_id, sequence, quality_scores, sequence_length
        @Claude
    """
    records = []
    
    for i, (header, sequence, plus_line, quality) in enumerate(parse_fastq_gz(filepath)):
        if max_records and i >= max_records:
            break
            
        # Extract sequence ID (remove @ prefix and any description)
        seq_id = header[1:].split()[0]
        
        records.append({
            'sequence_id': seq_id,
            'sequence': sequence,
            'quality_scores': quality,
            'sequence_length': len(sequence)
        })
    
    return pd.DataFrame(records)


def translate(dna_seq):
    '''function to translate DNA to protein sequence
    
    From: https://www.geeksforgeeks.org/dna-protein-python-3/'''
    
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(dna_seq)%3 == 0:
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i + 3]
            protein+= table[codon]
    else:
        protein = "frameshift"
    return protein

def safe_split(seq, front_flank, rear_flank):
        try: 
            front_split = seq.split(front_flank)
            if len(front_split) > 1:
                rear_split = front_split[1].split(rear_flank)
                if len(rear_split) > 1:
                    return rear_split[0]
                else:
                    return None
            else:
                return None
        except:
            return "lizard"
    
def safe_translate(seq):
    try:
        return translate(seq)
    except:
        return "lizard"

def main():
    parser = ArgumentParser(description="Process FASTQ.gz file and extract sequences.")
    parser.add_argument("fastq_gz_dir", type=str, help="Path to the directory containing (merged) FASTQ.gz files")
    parser.add_argument("outdir", type=str, help="Output directory for results")
    parser.add_argument("--save_seqs", action="store_true", dest="save_seqs", help="Whether to save all intermediate sequences to CSV")
    parser.add_argument('--designs', nargs='?', const=None, default=None, help="filter for seqs from ordered designs df passed as optional argument")
    args = parser.parse_args()
    fastq_gz_dir = args.fastq_gz_dir
    outdir = args.outdir
    save_seqs = args.save_seqs
    
    os.makedirs(outdir, exist_ok=True)

    print(f"Processing FASTQ.gz files in directory: {fastq_gz_dir}")
    print(f"Output directory: {outdir}")

    dfs_to_concat = []

    for fastq_file in glob.glob(os.path.join(fastq_gz_dir, "*.gz")):

        print(f"Processing file: {os.path.basename(fastq_file)}")
        file_id = os.path.basename(fastq_file).split(".")[0]
        
        try:
            # Read all sequences (remove max_records parameter to read entire file)
            df = fastq_gz_to_dataframe(fastq_file)
            
            print(f"Loaded {len(df)} sequences")
            # print("\nDataFrame info:")
            # print(df.info())
            
            # print("\nFirst few records:")
            # print(df.head())
            
            # print(f"\nSequence length statistics:")
            # print(df['sequence_length'].describe())
            
            # Optional: Save to CSV
            if save_seqs:
                df.to_csv(os.path.join(outdir, f"{file_id}_all_sequences.csv"), index=False)
            
        except FileNotFoundError:
            print(f"File {fastq_file} not found. Please check the file path.")
        except Exception as e:
            print(f"Error processing file: {e}")

        flank1 = 'tcggcttcgcatatg'
        flank2 = 'ctcgagggtggaggt'

        flank1 = flank1.upper()
        flank2 = flank2.upper()

        #save all seqs missing both flanks to separate df
        no_flank_df = df[~df['sequence'].str.contains(flank1) | ~df['sequence'].str.contains(flank2)].copy()
        no_flank_df.to_csv(os.path.join(outdir, f"{file_id}_sequences_missing_flanks.csv"), index=False)
        
        start = time.time()
        df['seq_to_translate'] = df['sequence'].apply(lambda x: safe_split(x, flank1, flank2))
        df['protein_sequence'] = df['seq_to_translate'].apply(lambda x: safe_translate(x) if x not in [None, "lizard"] else x)
        end = time.time()

        # print(f"Time to process {len(df)} sequences: {end - start:.2f} seconds, average {((end - start)/len(df)):.5f} seconds/sequence")
        
        seq_count_df = df['protein_sequence'].value_counts().reset_index()
        seq_count_df.columns = ['protein_sequence', 'count']
        # print(seq_count_df.head(20))
        if save_seqs:
            df.to_csv(os.path.join(outdir, f"{file_id}_all_sequences_with_protein.csv"), index=False)
        seq_count_df.to_csv(os.path.join(outdir, f"{file_id}_counts.csv"), index=False)

        #search for seqs which match ordered designs
        if args.designs is not None:
            designed_seqs_df = pd.read_csv(args.designs)
            designed_seqs_df = designed_seqs_df.rename(columns={'Name':'design'})

            #filter seq_count_df for seqs which match designs
            matching_seqs_df = seq_count_df[seq_count_df['protein_sequence'].isin(designed_seqs_df['Sequence'])]
            matching_seqs_df.to_csv(os.path.join(outdir, f"{file_id}_matching_sequences.csv"), index=False)

            dfs_to_concat.append(matching_seqs_df)
        
        print(f"Total reads: {df.shape[0]}")
        print(f"Number of seqs mising flanks: {no_flank_df.shape[0]}")
        print(f"Number of seqs with error in translation: {(df['protein_sequence'] == 'lizard').sum()}")
        print(f"Number of counts that map to library: {matching_seqs_df['count'].sum()}\n")
    
    if len(dfs_to_concat) != 0:
        print(f"Concatenating count tables for {len(dfs_to_concat)} files with matching sequences...")
        counts_df = pd.concat(dfs_to_concat, ignore_index=True, axis=1)
        counts_df.to_csv(os.path.join(outdir, f"count_table.csv"), index=False)
    else:
        print(f"Mapped sequences option not specified or no matching sequences found, skipping count table concatenation.\n")

# Example usage
if __name__ == "__main__":
    main()

    