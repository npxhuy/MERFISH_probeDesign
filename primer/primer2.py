import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from primer3 import calc_hairpin, calc_homodimer

def filter_and_trim_fasta(input_file, output_file):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            for i in range(len(record)-19):  # trim to 20 nt
                trimmed_seq = record.seq[i:i+20]

                # Check for hairpin formation and self-dimerization
                hairpin_tm = calc_hairpin(str(trimmed_seq)).tm
                homodimer_tm = calc_homodimer(str(trimmed_seq)).tm
                if hairpin_tm > 45 or homodimer_tm > 45:
                    continue

                # Remove sequences containing "CCCC" or "GGGG"
                if "CCCC" in trimmed_seq or "GGGG" in trimmed_seq:
                    continue

                # Remove sequences that do not have "C" or "G" in the last 5 characters
                if not ("C" in trimmed_seq[-5:] or "G" in trimmed_seq[-5:]):
                    continue

                # Only allow up to 3 'C' and 'G' bases in the last 5 bases
                if trimmed_seq[-5:].count('C') + trimmed_seq[-5:].count('G') > 3:
                    continue

                # Remove sequences with a melting temperature not between 59 and 60 degrees Celsius
                if not (59 <= mt.Tm_NN(trimmed_seq) <= 60):
                    continue

                # Remove sequences with a GC content less than 50% or greater than 60%
                gc_content = 100.0 * (trimmed_seq.count('G') + trimmed_seq.count('C')) / len(trimmed_seq)
                if gc_content < 35 or gc_content > 65:
                    continue

            

                # Write the trimmed sequence to the output file
                record.seq = trimmed_seq
                SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    input_file = sys.argv[1]  # get the input file name from the command-line arguments
    output_file = sys.argv[2]  # get the output file name from the command-line arguments
    filter_and_trim_fasta(input_file, output_file)