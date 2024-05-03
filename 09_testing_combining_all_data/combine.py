import os
import pandas as pd 

# Open the list of microbe names and read the contents
probe_path = '03_selected_probe'
name_df = pd.read_csv('name.txt', 
                   names=['specie'])

# Open readout seq file
readout_path = '02_readout_seq'
readout_df = pd.read_csv(os.path.join(readout_path, 
                                   'readout.txt'),
                                   sep = '\t', 
                                   names = ['position', 'sequence'])

# Open the codebook file
codebook_path = '01_codebook'
codebook_df = pd.read_csv(os.path.join(codebook_path, 
                                   'out.txt'),
                                   sep = '\t', 
                                   names = ['specie', 'positions','code'])

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
        

    '''
    pos1, pos2, pos3, pos4 = positions # Convert the list elements to integers

    # FROM POSITION, TO READOUT_SEQ
    readout1 = readout_df.loc[readout_df['position'] == pos1, 'sequence'].values[0]
    readout2 = readout_df.loc[readout_df['position'] == pos2, 'sequence'].values[0]
    readout3 = readout_df.loc[readout_df['position'] == pos3, 'sequence'].values[0]
    readout4 = readout_df.loc[readout_df['position'] == pos4, 'sequence'].values[0]

    
    
    #Create new columns for the readouts and initialize them with readout1,2,3,4 strings
    for i in range(1, 5):
        probe_df[f'Readout_{i}'] = locals()[f'readout{i}']
    
    # Remove each readout from the list of readout column
    for i in range(len(probe_df)):
        readout_column = (i % 4) + 1
        probe_df.loc[i, f'Readout_{readout_column}'] = ""
        '''
    
    # Create new columns for the readouts and initialize them with empty strings
    # All these conditions just to distribute the readout to the quarter of the probe_df equally
    # Can not find the better way to do this, yet.
    # Found out, the above for loop
    '''
    quarter_size = len(probe_df) // 4
    remainder = len(probe_df) % 4
    print (quarter_size, remainder)

    for i in range(1, 5):
        probe_df[f'Readout_{i}'] = ''

    if remainder == 1:
        probe_df.loc[quarter_size + 1 :, 'Readout_1'] = readout1

        probe_df.loc[:quarter_size, 'Readout_2'] = readout2
        probe_df.loc[quarter_size*2 + 1:, 'Readout_2'] = readout2

        probe_df.loc[:quarter_size*2 - 1 + 1, 'Readout_3'] = readout3
        probe_df.loc[quarter_size*3 + 1:, 'Readout_3'] = readout3

        probe_df.loc[:quarter_size*3 - 1 +1 , 'Readout_4'] = readout4
    
    elif remainder == 2:
        
        probe_df.loc[quarter_size + 1 :, 'Readout_1'] = readout1

        probe_df.loc[:quarter_size, 'Readout_2'] = readout2
        probe_df.loc[quarter_size*2 + 1 + 1:, 'Readout_2'] = readout2

        probe_df.loc[:quarter_size*2 - 1 + 1 + 1, 'Readout_3'] = readout3
        probe_df.loc[quarter_size*3 + 1 + 1:, 'Readout_3'] = readout3

        probe_df.loc[:quarter_size*3 - 1 + 1 + 1 , 'Readout_4'] = readout4
    
    elif remainder == 3:
    
        probe_df.loc[quarter_size + 1 :, 'Readout_1'] = readout1

        probe_df.loc[:quarter_size, 'Readout_2'] = readout2
        probe_df.loc[quarter_size*2 + 1 + 1:, 'Readout_2'] = readout2

        probe_df.loc[:quarter_size*2 - 1 + 1 + 1, 'Readout_3'] = readout3
        probe_df.loc[quarter_size*3 + 1 + 1 + 1:, 'Readout_3'] = readout3

        probe_df.loc[:quarter_size*3 - 1 + 1 + 1 + 1 , 'Readout_4'] = readout4

    else: 
        # Fill the 'Readout' columns with the readout values, leaving the corresponding quarter of each column empty
        probe_df.loc[quarter_size:, 'Readout_1'] = readout1

        probe_df.loc[:quarter_size - 1, 'Readout_2'] = readout2
        probe_df.loc[quarter_size*2:, 'Readout_2'] = readout2

        probe_df.loc[:quarter_size*2 - 1, 'Readout_3'] = readout3
        probe_df.loc[quarter_size*3 :, 'Readout_3'] = readout3

        probe_df.loc[:quarter_size*3 - 1 , 'Readout_4'] = readout4
        
    '''

    # Create new columns for the readouts and initialize them with empty strings
    # Calculate the size of each quarter

    # Create new columns for the readouts and initialize them with empty strings
    
    # Add 'Chr' and 'Seq' columns to combine_df
    combine_df = combine_df.append(probe_df.assign(Specie=specie, 
                                               Target=probe_df['Seq'],
                                               Primer_1="TBU",
                                               Primer_2="TBU",
                                               Chr=probe_df['Chr'],
                                               Readout_1=probe_df['Readout_1'],
                                               Readout_2=probe_df['Readout_2'],
                                               Readout_3=probe_df['Readout_3'],
                                               Readout_4=probe_df['Readout_4'])[['Specie', 'Primer_1', 'Primer_2', 'Chr', 'Target', 'Readout_1', 'Readout_2', 'Readout_3', 'Readout_4']], ignore_index=True)
# Print the new DataFrame
print(combine_df)