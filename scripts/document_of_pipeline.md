1. Download data for Paintshop Pipeline
Require Fasta (fna) and GTF to run the pipeline
Go to NCBI and search for the bacteria, for example _Prevotella nigrescens_ ([Prevotella nigrescens - NCBI - NLM (nih.gov)](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/28133/))
Click on download, RefSeq only, check on GTF and Fasta, optional to name the zip file (recommend to do so)
cd When extract the path is like: prevotella_nigrescens/ncbi_dataset/data/"bacterial ID"/"fna + gtf files"
   
   
2. Extract the data because the path file is a mess
So I make two folders, one is lungs_microbes, contains all the raw extracted file and everything in it
One is extract_data where i take only the fna and gtf from all the messy path

Run this inside the lungs_microbes to make the same folder name in the extract_data
`ls | while read folder; do mkdir ../extract_data/$folder
And this to copy the file and rename it (i actually mv instead of cp, it's stupid, with this file it's small so no problem, but with bigger file then should link the file to another folder instead and having a system to note all these thing down)
`ls | while read folder; do mkdir ../extract_data/$folder; cp $folder/*/*/*/*.fna ../extract_data/$folder/$folder.fna; cp $folder/*/*/*/*.gtf ../extract_data/$folder/$folder.gtf; done`

3. The GTF files has some # so i have to remove it, run it inside the extract_data
`ls | while read folder; do cd $folder; cat *.gtf | grep -v "#" > temp.gtf && mv temp.gtf $folder.gtf; cd ..; done`

5. Run the paintshop pipeline on the extract_data
NOTE: remember to specify the temperature and length of probes in the loop_paintshop.sh when making config file, or else it will use the default parameter

6. Make Blast DB
Make directory for containing blast database for lateral blastn run
Run this in the extract_data folder
`ls | while read folder; do mkdir ../../blast_db/$folder; done  `

Combine the fna file, exclude one microbe to make db
`ls | while read folder; do ls | grep -v $folder | while read microbe; do cat $microbe/$microbe.fna >> ../for_blast/$folder/$folder.fasta ; done`


`ls | while read folder; do makeblastdb -in $folder/$folder.fasta -dbtype nucl -out ../../blast_db/$folder/$folder; done`
7. Move the probes to another folder
  `ls | while read folder; do mkdir ../../probes/$folder; done`

Filter the file, and wrote it as a fasta file cause blastn want it fasta
  `ls | while read folder; do cat $folder/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | cut -f4 | awk '{print ">\n"$0}' > ../../probes/$folder/$folder.fasta; done`

8. Run blastn 
`ls | while read folder; do mkdir ../../blastn/$folder; done`

Run blastn on the db
`ls | while read folder; do blastn -db ../blast_db/$folder/$folder -query $folder/$folder.fasta -word_size 15 -ungapped > ../blastn/$folder/$folder.txt`


9. Filter uniq probes
Run quite slow
`ls | while read folder; do grep "Sbjct" -B 2 $folder/$folder.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - ../probes/$folder/$folder.fasta | grep -v ">" > ../uniq_probes/uniq_$folder.txt; done`
   