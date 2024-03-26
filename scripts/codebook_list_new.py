from itertools import combinations

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def generate_codewords(length):
    # Generate all possible combinations of indices for placing four 1s
    indices = list(range(length))
    codewords = []
    for comb in combinations(indices, 4):
        codeword = ['0'] * length
        prev_index = None
        for index in comb:
            # Ensure no more than three consecutive 1s
            if prev_index is not None and index == prev_index + 1:
                break
            codeword[index] = '1'
            prev_index = index
        else:  # This block executes if the loop does not break
            codeword_str = ''.join(codeword)
            # Check Hamming distance
            valid = True
            for existing_codeword in codewords:
                if hamming_distance(codeword_str, existing_codeword) < 4:
                    valid = False
                    break
            if valid:
                codewords.append(codeword_str)
    return codewords

# Example usage:
length = int(input("Enter the length of the codewords: "))

# Generate codewords
codewords = generate_codewords(length)

print("Generated Codewords:")
for cw in codewords:
    print(cw)
print("Total number of codewords:", len(codewords))
