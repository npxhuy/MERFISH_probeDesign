# Read the file 'uniq_species_t10.txt'
cat uniq_species_t10.txt | 

# For each line in the file (each 'microbe')
while read -r microbe; do 

  # Replace spaces in the microbe name with underscores and assign to 'folder'
  folder=$(echo $microbe | tr ' ' '_'); 

  # Use the 'datasets' tool to download genome data for the microbe
  # The downloaded data is:
  # - Annotated
  # - Complete assembly level
  # - From the RefSeq source
  # - The latest assembly version
  # - Includes the genome and GTF files
  # The data is saved in a zip file named after the 'folder' variable
  datasets download genome taxon "$microbe" --annotated --assembly-level complete --assembly-source refseq --assembly-version latest --reference --include genome,gtf --filename $folder.zip; 

# End of the while loop
done