import os
import pandas as pd 

# Open the list of microbe names and read the contents
probe_path = '03_selected_probes'
#name_df = pd.read_csv('final_microbes.txt', 
#                   names=['specie'])

# Open readout seq file
readout_path = '02_readout_seq'
readout_df = pd.read_csv(os.path.join(readout_path, 
                                   'readout.txt'),
                                   sep = '\t', 
                                   names = ['position', 'sequence'])

# Open the codebook file
codebook_path = '01_codebook'
codebook_df = pd.read_csv(os.path.join(codebook_path, 
                                   'microbes_codebook.txt'),
                                   sep = '\t', 
                                   names = ['specie', 'positions','code'])

name_df = pd.read_csv(os.path.join(codebook_path, 
                                   'final_microbes.txt'), 
                                    names=['specie'])

# Create a new DataFrame with the specified column names
combine_df = pd.DataFrame(columns=['Specie', 'Chr', 'Primer_1', 
                                   'Readout_1', 'Readout_2', 'Target',
                                   'Readout_3', 'Readout_4', 'Primer_2'])

for index, row in name_df.iterrows():
    specie = row['specie']
 
    column_names = ['Chr', 'Start', 'Stop', 
                    'Seq', 'Tm', 'on-targ', 
                    'off-targ', 'repeat', 'max k-m', 
                    'strand', 'GeneID', 'Transcript']
    probe_df = pd.read_csv(os.path.join(probe_path, f'{specie}.txt'),
                       sep='\s+', #not \t because the separator of geneid and transcript is sth idk
                       names=column_names)
      
    # Fill the 'Readout' columns with the readout values, leaving the corresponding quarter of each column empty

    # POSITION FROM CODEBOOK
    positions_str = codebook_df.loc[codebook_df['specie'] == specie, 'positions'].values[0]
    positions = positions_str.strip('[]').split() # Remove the brackets and split the string into a list
    positions = [int(pos) for pos in positions] # Convert the list elements to integers
    

    # FROM POSITION, TO READOUT_SEQ
    readouts = [readout_df.loc[readout_df['position'] == pos, 'sequence'].values[0] for pos in positions]

    # Create new columns for the readouts and initialize them with readout1,2,3,4 strings
    for i in range(len(readouts)):
        probe_df[f'Readout_{i+1}'] = readouts[i]

    # Remove each readout from the list of readout column
    for i in range(len(probe_df)):
        readout_column = (i % len(readouts)) + 1
        probe_df.loc[i, f'Readout_{readout_column}'] = ""
        
    # Create new columns for the readouts and initialize them with empty strings
    # Calculate the size of each quarter

    # Create new columns for the readouts and initialize them with empty strings
    
    # Add 'Chr' and 'Seq' columns to combine_df
    ### PLEASE READ THIS COMMENT
    ### On pandas version 1, the method is called .append instead of ._append
    ### If you are using pandas version 2, then it's fine.
    ### If you are using pandas version 1, change it to combine_df.append
    # Concatenate the columns into a new 'Probe' column
    combine_df = combine_df._append(probe_df.assign(Specie=specie, 
                                               Target=probe_df['Seq'],
                                               Primer_1="CGTGTTAGTGGCCCGGGTCT",
                                               Primer_2="GGCCGCGACTAGGTAAGCCT",
                                               Chr=probe_df['Chr'],
                                               Readout_1=probe_df['Readout_1'],
                                               Readout_2=probe_df['Readout_2'],
                                               Readout_3=probe_df['Readout_3'],
                                               Readout_4=probe_df['Readout_4'])[['Specie', 'Primer_1', 'Primer_2', 'Chr', 'Target', 'Readout_1', 'Readout_2', 'Readout_3', 'Readout_4']], ignore_index=True)
# Print the new DataFrame
combine_df.to_csv('probe_table.tsv', sep = '\t', index=False)