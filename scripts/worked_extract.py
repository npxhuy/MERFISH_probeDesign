import sys
import pandas as pd

# Get the file paths from the command line arguments
gtf_path = sys.argv[1]  # Path to the GTF file
probe_location_path = sys.argv[2]  # Path to the probe location file
output_file_path = sys.argv[3]  # Path to the output file

# Load the data from the input files into pandas DataFrames
gtf = pd.read_csv(gtf_path, sep='\t', header=None, names=['ID', 'Start', 'End', 'Gene'])
probe_location = pd.read_csv(probe_location_path, sep='\t', header=None, names=['ID', 'Start', 'End'])

# Initialize a new column in probe_location with 'No Match' as the default value
probe_location['Gene'] = 'No Match'

# Iterate over the rows in gtf
for i, row in gtf.iterrows():
    # For each row in gtf, create a mask that is True for rows in probe_location where:
    # - The 'ID' matches the 'ID' in the current row of gtf
    # - The 'Start' is greater than or equal to the 'Start' in the current row of gtf
    # - The 'End' is less than or equal to the 'End' in the current row of gtf
    mask = (probe_location['ID'] == row['ID']) & (probe_location['Start'] >= row['Start']) & (probe_location['End'] <= row['End'])
    
    # Update the 'Gene' column in probe_location for the rows where the mask is True
    # Set the 'Gene' value to the 'Gene' value from the current row of gtf
    probe_location.loc[mask, 'Gene'] = row['Gene']

# Save the result to the output file
# The output file is a tab-separated values (TSV) file without an index column
probe_location.to_csv(output_file_path, sep='\t', index=False)