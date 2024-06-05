# List the contents of the current directory
ls | while read folder; do 

    # For each item, create a corresponding directory in ../05_blast_plus/01_fasta_combined/
    mkdir ../05_blast_plus/01_fasta_combined/$folder

    # List the contents of the current directory again
    ls | 

    # Exclude the current item
    grep -v $folder | 

    while read microbe; do 

        # Concatenate the .fna file from each other directory into a single file in the new directory
        cat $microbe/$microbe.fna >> ../05_blast_plus/01_fasta_combined/$folder/$folder.fasta

    done

done