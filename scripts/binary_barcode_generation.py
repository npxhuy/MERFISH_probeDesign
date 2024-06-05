import re
import sys

def extract_barcode_positions_and_length(html_file):
    barcode_positions = []
    barcode_length = 0

    with open(html_file, 'r') as file:
        for i, line in enumerate(file, start=1):
            if i == 5:
                # Extract the barcode length from line 5
                match = re.search(r'>C\((\d+),', line)
                if match:
                    barcode_length = int(match.group(1))
            elif i >= 6 and not line.startswith('</pre>'):
                barcode_positions.append(line.strip())
            elif line.startswith('</pre>'):
                break

    return barcode_positions, barcode_length

def generate_barcodes(input_file, output_file):
    # Extract the barcode positions and length from the HTML file
    barcode_positions, barcode_length = extract_barcode_positions_and_length(input_file)

    with open(output_file, "w") as file:
        for line in barcode_positions:
            # Initialize the barcode as a list of zeros
            barcode = [0] * barcode_length

            # Split the line into numbers
            positions = list(map(int, line.split()))

            # Set the positions in the barcode to 1
            for pos in positions:
                barcode[pos - 1] = 1  # Subtract 1 because positions are 1-based

            # Convert the barcode to a string and write it to the output file
            file.write(''.join(map(str, barcode)) + "\n")

# Call the function with the names of your input and output files
generate_barcodes(sys.argv[1], sys.argv[2])