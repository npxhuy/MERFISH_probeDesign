from itertools import permutations
import multiprocessing
import time

# Define the sequence
sequence = 'G'*5 + 'A'*3 + 'T'*3

# Define a function to generate and filter permutations
def generate_and_filter_permutations(start_char):
    # Generate all permutations of the sequence that start with the specified character
    perms = set(p for p in permutations(sequence) if p[0] == start_char)
    # Filter out permutations that contain four or more consecutive 'G's
    filtered_perms = [p for p in perms if 'GGGG' not in ''.join(p)]
    # Return the filtered permutations
    return filtered_perms

# Record the start time
start_time = time.time()

# Create a pool of worker processes
with multiprocessing.Pool() as pool:
    # Run the three tasks in parallel
    tasks = ['A', 'T', 'G']
    results = pool.map(generate_and_filter_permutations, tasks)

# Open the output file
with open('ALL_readout_output.txt', 'w') as f:
    # Write the permutations to the file
    for result in results:
        for p in result:
            # Convert the permutation to a string
            s = ''.join(p)
            f.write(s + '\n')

# Record the end time
end_time = time.time()

# Calculate the elapsed time
elapsed_time = end_time - start_time
print(f'Elapsed time: {elapsed_time} seconds')