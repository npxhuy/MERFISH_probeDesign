# Change the csv probes file to fasta format
awk '{print ">\n"$0}' filename > newfile

# Bowtie2 pipeline

# Create bowtie2 index
bowtie2-build pn.fna bowtie2/pn

# Align reads to reference genome
# end to end alignment as default mode
bowtie2 -f --local -D 15 -R 2 -N 1 -L 15 -i S,1,0.75 -x bowtie2/pn -U sa.fasta -S sa_pn_3.sam

# Explanation of bowtie2 flag and terms
# Seed: short, fixed-length sequence used to find alignment
# -L indicates the seed length
# e.g our reads (probes) are 30bp long, and we want the minimum discarding length to be 15bp
# from the 15bp seed, it starts to find the alignment
# -D indicates the maximum number of seed extension attempts
# -N indicates the maximum number of mismatches in a seed alignment
# -R indicates the maximum number of re-seed attempts
# to explain like I am five, it means: there are too many possible alignment and we dont want to waste time do it over again
# let's try to look for sth else now


# Take only the mapped reads
# Replace != with == for unmapped reads -> can filter the unmapped reads for later
awk -F'\t' '$3 != "*"' sa_pn_3.sam


## BLASTN pipeline
# could it be better?

# Create blast database


blastn -db blast/pn -query sa.fasta -word_size 15 -ungapped > sa_pn_blastn_2.txt
#both gapped and ungapped result in the same result, but prob will go for ungapped

# Extract the sequence that matched with the reference genome
# And filter the sequence that matched with the reference genome against the original fasta file
# to have the unique non-matching sequence
grep "Sbjct" -B 2 sa_pn_blastn_2.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - sa.txt > uniq_sa.txt

grep "Sbjct" -B 2 prevotella_bryantii.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - prevotella_bryantii.fasta | grep -v ">" > uniq_prevotella_bryantii.txt


ls | while read folder; do grep "Sbjct" -B 2 $folder/$folder.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - ../probes/$folder/$folder.fasta | grep -v ">" > ../uniq_probes/uniq_$folder.txt; done

ls | while read folder; do cat $folder/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | cut -f4 | awk '{print ">\n"$0}' > ../../probes/$folder/$folder.fasta



ls | while read folder; do blastn -db ../blast_db/$folder/$folder -query $folder/$folder.fasta -word_size 15 -ungapped > ../blastn/$folder/$folder.txt

