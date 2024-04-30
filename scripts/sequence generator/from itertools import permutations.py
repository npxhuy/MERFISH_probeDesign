from itertools import permutations
import multiprocessing as mp
import time
# Define the sequence
sequence = 'G'*5 + 'A'*4 + 'T'*3


# Define a function to generate and filter permutations
def generate_and_filter_permutations(start_char):
    # Generate all permutations of the sequence that start with the specified character
    perms = set(p for p in permutations(sequence) if p[0] == start_char)
    # Filter out permutations that contain four or more consecutive 'G's
    filtered_perms = [p for p in perms if 'GGGG' not in ''.join(p)]
    # Return the filtered permutations
    return filtered_perms

if __name__ == '__main__':
    
    pA = mp.Process(target=generate_and_filter_permutations, args=('A',))
    pT = mp.Process(target=generate_and_filter_permutations, args=('T',))
    pG = mp.Process(target=generate_and_filter_permutations, args=('G',))
    
    start_time = time.time()
    pA.start()
    pT.start()
    pG.start()
    end_time = time.time()


    # Calculate the elapsed time
    elapsed_time = end_time - start_time
    print(f'Elapsed time: {elapsed_time} seconds')