# List the contents of the current directory
ls | while read folder; do 

    # For each item, create a corresponding directory
    mkdir /home/npxhuy/03_data/05_blast/02_makedb_01/$folder

    # Run the makeblastdb command on the .fasta file in each directory
    # The output database is stored in the new directory
    /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/makeblastdb \
        -in $folder/$folder.fasta \
        -dbtype nucl \
        -out /home/npxhuy/03_data/05_blast/02_makedb_01/$folder/$folder

done