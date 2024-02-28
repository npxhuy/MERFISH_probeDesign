import sys
import pandas as pd

# Accept the file paths as command line arguments
gtf_path = sys.argv[1]
probe_location_path = sys.argv[2]
output_file_path = sys.argv[3]  # New line

# Open the files
gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None)
location = pd.read_csv(probe_location_path, sep='\t', header=None)

# Open the output file
with open(output_file_path, 'w') as output_file:  # Modified line
    j=0
    for i in range(location.shape[0]):
        while j < gtf.shape[0]:
            if gtf.iloc[j, 0] == location.iloc[i, 0]:
                
                if gtf.iloc[j, 1] > location.iloc[i, 1] and gtf.iloc[j, 2] > location.iloc[i, 2]:
                    output_file.write('No Match' + '\n')
                    break

                elif gtf.iloc[j, 1] <= location.iloc[i, 1] and gtf.iloc[j, 2] >= location.iloc[i, 2]:
                    # Write the fourth column of the current row in gtf to the output file
                    output_file.write(str(gtf.iloc[j, 3])'\n')
                    break

                elif gtf.iloc[j, 1] < location.iloc[i, 1] and gtf.iloc[j, 2] < location.iloc[i, 2]:
                    if j+1 < gtf.shape[0] and gtf.iloc[j+1 , 1] == location.iloc[i+1, 1]:
                        j += 1
                        break
                    elif j+1 < gtf.shape[0] and gtf.iloc[j+1 , 1] != location.iloc[i+1, 1]:
                        j += 1

            
            
                elif gtf.iloc[j,0] != location.iloc[i, 0]:
                    if j+1 < gtf.shape[0] and i+1 < location.shape[0] and gtf.iloc[j+1,0] != location.iloc[i+1, 0]:
                        j += 1
                    elif j+1 < gtf.shape[0] and i+1 < location.shape[0] and gtf.iloc[j+1,0] == location.iloc[i+1, 0]:
                        break

            else:
                j += 1
