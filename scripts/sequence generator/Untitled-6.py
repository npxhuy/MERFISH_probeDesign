#!/bin/python
#SBATCH -A lu2023-7-68
#SBATCH -p gpua100i
#SBATCH -t 10-00:00:00
#SBATCH -N 3
#SBATCH--ntasks-per-node=8
#SBATCH --mail-user=ph4342ng-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH -J para
#SBATCH -o /home/npxhuy/lu2023-17-27/hy/thesis/03_data/07_readout/04_testing_code_speed/testing_new/para.out
#SBATCH -e /home/npxhuy/lu2023-17-27/hy/thesis/03_data/07_readout/04_testing_code_speed/testing_new/para.out
#SBATCH -D /home/npxhuy/lu2023-17-27/hy/thesis/03_data/07_readout/04_testing_code_speed/testing_new

from itertools import permutations
import multiprocessing as mp
import time

# Define the sequence
sequence = 'G'*10 + 'A'*5 + 'T'*5

# Define a function to generate and filter permutations
def generate_and_filter_permutations(start_char):
    # Generate all permutations of the sequence that start with the specified character
    perms = set(p for p in permutations(sequence) if p[0] == start_char)
    # Filter out permutations that contain four or more consecutive 'G's
    filtered_perms = [p for p in perms if 'GGGG' not in ''.join(p)]
    # Write the filtered permutations to a file
    with open('Do_not_delete_output_10G_5A_5T.txt', 'a') as f:
        for p in filtered_perms:
            f.write(''.join(p) + '\n')

if __name__ == '__main__':
    
    pA = mp.Process(target=generate_and_filter_permutations, args=('A',))
    pT = mp.Process(target=generate_and_filter_permutations, args=('T',))
    pG = mp.Process(target=generate_and_filter_permutations, args=('G',))
    
    print('Multiprocessing...')
    start_time = time.time()
    pA.start()
    pT.start()
    pG.start()
    pA.join()
    pT.join()
    pG.join()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'Elapsed time: {elapsed_time} seconds')
    print('Done!')