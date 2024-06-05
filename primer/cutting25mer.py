import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def cut_sequences():
    input_file = sys.argv[1]  # First command-line argument
    output_file = sys.argv[2]  # Second command-line argument

    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            # Check if the sequence length is at least 25
            if len(record.seq) >= 25:
                # Cut the sequence into overlapping sequences of length 20
                for i in range(len(record.seq) - 19):
                    new_seq = record.seq[i:i+20]
                    if len(new_seq) == 20:
                        new_record = SeqRecord(new_seq, id=record.id, description="")
                        SeqIO.write(new_record, output_handle, "fasta")

# Call the function
cut_sequences()