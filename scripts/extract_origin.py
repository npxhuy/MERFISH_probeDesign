"""
This script reads two tab-separated files and compares their rows based on certain conditions.

The script accepts two command line arguments:
1. The path to the gtf (gtf_path)
2. The path to the second file contains the extracted location of the probes (probe_location_path)

The files are read into pandas DataFrames. The script then iterates over each row in the second DataFrame (location). For each row in the second DataFrame, it iterates over each row in the first DataFrame (gtf) until it finds a row that meets one of the following conditions:

1. If the first column of the current row in the first DataFrame matches the first column of the current row in the second DataFrame, and the third column of the current row in the first DataFrame is less than or equal to the third column of the current row in the second DataFrame, it moves on to the next row in the first DataFrame.

2. If the first column of the current row in the first DataFrame matches the first column of the current row in the second DataFrame, and the second column of the current row in the second DataFrame is greater than or equal to the second column of the current row in the first DataFrame, and the third column of the current row in the second DataFrame is less than or equal to the third column of the current row in the first DataFrame, it prints the fourth column of the current row in the first DataFrame and breaks the loop.

If neither condition is met, it moves on to the next row in the first DataFrame.

The script continues until it has iterated over all rows in the second DataFrame.
"""

import pandas as pd
import sys
from datetime import datetime

# Accept the file paths as command line arguments
gtf_path = sys.argv[1]
probe_location_path = sys.argv[2]

# Open the files
gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None)
location = pd.read_csv(probe_location_path, sep='\t', header=None)

# Record the start time
start_time = datetime.now()

# Open the output file
with open('output2.txt', 'w') as output_file:
    for i in range(location.shape[0]):
        j = 0
        while j < gtf.shape[0]:
            if gtf.iloc[j, 0] == location.iloc[i, 0]:
                if gtf.iloc[j, 2] <= location.iloc[i, 2]:
                    j += 1
                elif gtf.iloc[j, 1] <= location.iloc[i, 1] and gtf.iloc[j, 2] >= location.iloc[i, 2]:
                    # Write the fourth column of the current row in gtf to the output file
                    output_file.write(str(gtf.iloc[j, 3]) + '\n')
                    break
                elif gtf.iloc[j, 1] >= location.iloc[i, 2] and gtf.iloc[j, 2] >= location.iloc[i, 2]:
                    output_file.write('no_info' + '\n')
                    break 
            else:
                j += 1

# Print the elapsed time
print('Elapsed time: ', datetime.now() - start_time)