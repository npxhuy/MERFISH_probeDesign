def generate_barcodes(input_file, output_file):
    # Ask the user for the barcode length
    barcode_length = int(input("Please enter the barcode length: "))

    with open(input_file, "r") as file:
        lines = file.readlines()

    with open(output_file, "w") as file:
        for line in lines:
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
generate_barcodes("input.txt", "output.txt")