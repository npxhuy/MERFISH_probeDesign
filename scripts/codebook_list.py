from itertools import combinations

def hamming_distance(codeword1, codeword2):
    return sum(c1 != c2 for c1, c2 in zip(codeword1, codeword2))

def generate_codewords(length, min_hamming_distance):
    # Generate all possible combinations of indices for placing the 1s
    indices = list(range(length))
    codewords = []
    for comb in combinations(indices, 4):
        codeword = ['0'] * length
        for index in comb:
            codeword[index] = '1'
        codeword_str = ''.join(codeword)
        # Check if the generated codeword meets the Hamming distance criteria
        valid = True
        for existing_codeword in codewords:
            if hamming_distance(codeword_str, existing_codeword) < min_hamming_distance:
                valid = False
                break
        if valid:
            codewords.append(codeword_str)
    return codewords

# Example usage:
length = int(input("Enter the length of the codewords: "))
min_hamming_distance = int(input("Enter the minimum Hamming distance: "))

# Generate codewords
codewords = generate_codewords(length, min_hamming_distance)

filename = '{}bit_HD{}_{}.txt'.format(length, min_hamming_distance, len(codewords))

# Open the file in write mode
with open(filename, 'w') as f:
    # Write the codewords to the file
    for cw in codewords:
        f.write(cw + '\n')