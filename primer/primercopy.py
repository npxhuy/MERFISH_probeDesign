import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from primer3 import calc_hairpin, calc_homodimer
import random

def filter_and_trim_fasta(input_file, output_file):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            # Remove sequences containing "CCCC" or "GGGG"
            
            if "CCCC" in record.seq or "GGGG" in record.seq:
                continue

            # Remove sequences that do not have "C" or "G" in the last 5 characters
            if not ("C" in record.seq[-5:] or "G" in record.seq[-5:]):
                continue

            # Only allow up to 3 'C' and 'G' bases in the last 5 bases
            if record.seq[-5:].count('C') + record.seq[-5:].count('G') > 3:
                continue

            # Remove sequences with a melting temperature not between 58 and 60 degrees Celsius
            if not (58 <= mt.Tm_NN(record.seq) <= 60):
                continue

            # Remove sequences with a GC content less than 50% or greater than 60%
            gc_content = 100.0 * (record.seq.count('G') + record.seq.count('C')) / len(record.seq)
            if gc_content < 35 or gc_content > 65:
                continue

            # Check for hairpin formation and self-dimerization
            hairpin_tm = calc_hairpin(str(record.seq)).tm
            homodimer_tm = calc_homodimer(str(record.seq)).tm
            if hairpin_tm < 40 or homodimer_tm < 40:
                continue
            
            

            SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    input_file = sys.argv[1]  # get the input file name from the command-line arguments
    output_file = sys.argv[2]  # get the output file name from the command-line arguments
    filter_and_trim_fasta(input_file, output_file)