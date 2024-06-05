from primer3 import calc_hairpin, calc_homodimer
seq = "GGCCGCGACTAGGTAAGCCT"

print(calc_hairpin(seq, max_loop=10))
print(calc_homodimer(seq, max_loop=10))

hairpin_tm = calc_hairpin(seq).tm
homodimer_tm = calc_homodimer(seq).tm
if hairpin_tm < 50 and homodimer_tm < 50:
    print("Yes lower than 50")
else:
    print("No greater than 50")
                    