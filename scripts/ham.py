def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

with open("a.txt") as f1, open("b.txt") as f2, open("pb_filter.txt", "w") as f3, open("pn_filter.txt", "w") as f4:
    list1 = [line.strip() for line in f1]
    list2 = [line.strip() for line in f2]

    # Filter list1
    filtered_list1 = [elem1 for elem1 in list1 if all(hamming_distance(elem1, elem2) > 15 for elem2 in list2)]
    for elem1 in filtered_list1:
        f3.write(elem1 + "\n")

    # Filter list2
    filtered_list2 = [elem2 for elem2 in list2 if all(hamming_distance(elem2, elem1) > 15 for elem1 in list1)]
    for elem2 in filtered_list2:
        f4.write(elem2 + "\n")