# REMEMBER TO DOWNLOAD BIOPYTHON
# pip3 install biopython



# In case of SSL error, try to install certificate.
# https://gist.github.com/marschhuynh/31c9375fc34a3e20c2d3b9eb8131d8f3#file-install-certificates-command

# Load biopython module.
from Bio import Entrez 
from Bio import SeqIO

# Modify line 18 to your own email address.
# Modify line 21 to your own search term.

# Set the email address to use for NCBI requests.
def download_seq(species_list):
    Entrez.email = "ph4342ng-s@student.lu.se"

    for species in species_list: # Loop through the species list.
        search_term = f"{species} AND 16s ribosomal RNA gene AND 1:1600[SLEN], NOT shotgun, NOT partial" # Search term for NCBI.
        handle = Entrez.esearch(db="nuccore", term=search_term, retmax=5) # Search NCBI.
        record = Entrez.read(handle) # Read the search results.
        handle.close() # Close the handle.

        if int(record["Count"]) == 0: # If no results are found, print an error message.
            print(f"No complete genome found for {species}") # Print error message.
        else:    
            seq_id = record["IdList"] # Get the sequence ID.
            #Fasta file
            handle_fasta = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text") # Get the fasta file.
            sequences = SeqIO.parse(handle_fasta, "fasta") # Parse the fasta file.
            sequence_list = list(sequences) # Convert the fasta file to a list.
            output_fasta = f"{species}.fasta" # Set the output file name.
            SeqIO.write(sequence_list, output_fasta, "fasta") # Write the fasta file.
            handle_fasta.close() # Close the handle.
            num_sequences = len(sequence_list) # Get the number of sequences.
            print(f"Downloaded {num_sequences} sequences for {species}") # Print the number of sequences downloaded.

file_path = 'species_list.txt' # Path to the species list.

with open(file_path,'r') as file: # Open the species list.
    species_list = [line.strip() for line in file] # Read the species list.

download_seq(species_list) # Run the function to download the sequences.