# This script is used to match and merge data from two tab-separated files based on specific conditions.
# It uses the pandas library to handle the data.

# Import necessary libraries: sys for command-line arguments and pandas for data manipulation.

# Retrieve the file paths from the command-line arguments. The first two arguments are the input files and the third argument is the output file.

# Load the data from the input files into pandas DataFrames (df1 - the gtf file and df2 - the probe location file). The data is assumed to be tab-separated. The column names are specified for each DataFrame.

# Initialize a new column in df2 named 'GeneID' and 'Transcript', and fill it with the value 'No_Match'. This column will be updated later with matching values from df1.

# Iterate over each row in df1. For each row, find rows in df2 where 'Chr' matches, 'Start' is greater than or equal to the 'Start' in df1, and 'Stop' is less than or equal to the 'Stop' in df1.

# For each matching row in df2, update the 'GeneID' and 'Transcript' columns with the corresponding values from df1.

# After all rows in df1 have been processed, save df2 to a new tab-separated file. The file path is specified by the third command-line argument.

# This script is typically used for processing genomic data, where df1 contains gene annotations (chromosome, start and stop positions, gene ID, and transcript) and df2 contains some genomic features or measurements (chromosome, start and stop positions, sequence, melting temperature, on-target score, off-target score, repeat score, maximum k-mer score, and strand). The script matches the features in df2 to the genes in df1 based on their positions.

import sys
import pandas as pd

# Get the file paths from the command line arguments
file1_path = sys.argv[1]
file2_path = sys.argv[2]
output_file_path = sys.argv[3]

# Load the data
df1 = pd.read_csv(file1_path, sep='\t', header=None, names=['Chr', 'Start', 'Stop', 'GeneID', 'Transcript'])
df2 = pd.read_csv(file2_path, sep='\t', index_col=False, header=None, names=['Chr', 'Start', 'Stop', 'Seq', 'Tm', 'on-targ', 'off-targ', 'repeat', 'max k-m', 'strand'])

# Initialize a new column in df2
df2['GeneID'] = 'No_Match'
df2['Transcript'] = 'No_Match'

# Iterate over df1
for i, row in df1.iterrows():
    # Find matching rows in df2
    mask = (df2['Chr'] == row['Chr']) & (df2['Start'] >= row['Start']) & (df2['Stop'] <= row['Stop'])
    # Update the 'Gene' column in df2 for matching rows
    df2.loc[mask, 'GeneID'] = row['GeneID']
    df2.loc[mask, 'Transcript'] = row['Transcript']
        
# Save the result
df2.to_csv(output_file_path, sep='\t', index=False)