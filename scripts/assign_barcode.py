import sys
import random

# Get the input, output, and names file names from the command line arguments
barcode_filename = sys.argv[1]
microbes_filename = sys.argv[2]
output_filename = sys.argv[3]

# Read the names from the names file
with open(microbes_filename, 'r') as microbes_file:
    names = [line.strip() for line in microbes_file]

# Create a list to store the binary strings and positions
binary_strings_and_positions = []

# Open the input file
with open(barcode_filename, 'r') as barcode_file:
    # Read each line from the input file
    for line in barcode_file:
        # Remove the newline character
        binary_string = line.strip()
        # Get the positions of the '1's
        positions = [i+1 for i, bit in enumerate(binary_string) if bit == '1']
        # Convert the list of positions to a string and remove the commas
        positions_str = str(positions).replace(',', '')
        # Add the binary string and positions to the list
        binary_strings_and_positions.append((binary_string, positions_str))

# Randomly shuffle the binary strings and positions
random.shuffle(binary_strings_and_positions)

# Open the output file
with open(output_filename, 'w') as output_file:
    # Assign a name, binary string, and positions to each line in the output file
    for i in range(max(len(names), len(binary_strings_and_positions))):
        # Get the name for this line, or "NO_ASSIGNED" if there are no more names
        name = names[i] if i < len(names) else "NO_ASSIGNED"
        # Get the binary string and positions for this line, or "NO_BINARY_STRING" and "NO_POSITIONS" if there are no more binary strings and positions
        binary_string, positions_str = binary_strings_and_positions[i] if i < len(binary_strings_and_positions) else ("NO_BINARY_STRING", "NO_POSITIONS")
        # Write the name, binary string, and positions to the output file
        output_file.write(f'{name}\t{positions_str}\t{binary_string}\n')