import sys
import pandas as pd

# Get the file paths from the command line arguments
file1_path = sys.argv[1] # The probes' location file
file2_path = sys.argv[2] # The readouts' location file
output_file_path = sys.argv[3]

# Load the data
df1 = pd.read_csv(file1_path, sep=' ', header=None, names=['ID', 'Start', 'End'])
df2 = pd.read_csv(file2_path, sep=' ', header=None, names=['ID', 'Start', 'End'])

# Initialize a new column in df1
df1['Match'] = 'no'

# Iterate over df2
for i, row in df2.iterrows():
    # Find matching rows in df1
    mask = (df1['ID'] == row['ID']) & (df1['Start'] <= row['Start']) & (df1['End'] >= row['End'])
    # Update the 'Match' column in df1 for matching rows
    df1.loc[mask, 'Match'] = 'yes'

# Save the result
df1.to_csv(output_file_path, sep=' ', index=False)