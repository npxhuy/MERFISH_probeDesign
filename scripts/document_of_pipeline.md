1. Download data for Paintshop Pipeline

Require Fasta (fna) and GTF to run the pipeline
Go to NCBI and search for the bacteria, for example _Prevotella nigrescens_ ([Prevotella nigrescens - NCBI - NLM (nih.gov)](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/28133/))
Click on download, RefSeq only, check on GTF and Fasta, optional to name the zip file (recommend to do so)
cd When extract the path is like: prevotella_nigrescens/ncbi_dataset/data/"bacterial ID"/"fna + gtf files"

Nevermind: i found a way to download them quickly with the datasets tool from NCBI:

This is a list of microbes name and their taxon id

`cat microbes_list.tsv | while IFS=$'\t' read -r FIRST SECOND; do id=$(echo $SECOND | cut -d "/" -f 6); name=$(echo $FIRST | tr '[:upper:] ' '[:lower:]_'); echo $name $id; done`

This code to download them, replace the number after taxon with their id andm after file name, ofc made it in a loop.

`datasets download genome taxon 1857645 --assembly-source refseq --assembly-version latest --reference --include genome,gtf --filename Actinomyces_vulturis.zip`  

Combine two code, some species doesnt have complete genome on NCBI, optional to remove the -assembly-level

`cat microbes_list.tsv | while IFS=$'\t' read -r FIRST SECOND; do id=$(echo $SECOND | cut -d "/" -f 6); name=$(echo $FIRST | tr '[:upper:] ' '[:lower:]_'); datasets download genome taxon $id --annotated --assembly-level complete --assembly-source refseq --assembly-version latest --reference --include genome,gtf --filename $name.zip; done`

*Known issue* with the following three species, when run in the loop the resulted in this error
Error: Get "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon_suggest/2479840?tax_rank_filter=higher_taxon&taxon_resource_filter=TAXON_RESOURCE_FILTER_GENOME": read tcp 192.168.0.2:49951->130.14.29.110:443: read: connection reset by peer

But when run manually, it's normal, idk why\
Rothia koreensis \
Prevotella marseillensis\
Rothia aerolata


2. Extract the data because the path file is a mess.

So I make two folders, one is lungs_microbes, contains all the raw extracted file and everything in it
One is extract_data where i take only the fna and gtf from all the messy path

Run this inside the lungs_microbes to make the same folder name in the extract_data to copy the file and rename it (i actually mv instead of cp, it's stupid, with this file it's small so no problem, but with bigger file then should link the file to another folder instead and having a system to note all these thing down)

`ls | while read folder; do mkdir ../extract_data/$folder; cp $folder/*/*/*/*.fna ../extract_data/$folder/$folder.fna; cp $folder/*/*/*/*.gtf ../extract_data/$folder/$folder.gtf; done`

3. The GTF files has some # so i have to remove it, run it inside the extract_data

`for folder in */; do for file in "$folder"*.gtf; do grep -v "#" "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"; done; done`



4. Run the paintshop pipeline on the extract_data

NOTE: remember to specify the temperature and length of probes in the loop_paintshop.sh when making config file, or else it will use the default parameter

Known issue: problem on running paintshop on the following species:\
alcanivorax_profundi\
campylobacter_troglodytis\
fusobacterium_perfoetens\
haemophilus_seminalis\
methylobacterium_isbiliense\
parvimonas_parva\
peptostreptococcus_stomatis\
streptococcus_phocae\

Run seperately and 4 among them produced the probes but I can tell sth it's off, note it here:\
methylobacterium_isbiliense\
streptococcus_phocae\
alcanivorax_profundi\
campylobacter_troglodytis\


cat alcanivorax_profundi/pipeline_output/03*/01*/* > alcanivorax_profundi_Balance.tsv

6. Make Blast DB
Make directory for containing blast database for lateral blastn run
Run this in the extract_data folder
`ls | while read folder; do mkdir ../../blast_db/$folder; done  `

Combine the fna file, exclude one microbe to make db
`ls | while read folder; do ls | grep -v $folder | while read microbe; do cat $microbe/$microbe.fna >> ../blast_fna/$folder/$folder.fasta ; done; done`


`ls | while read folder; do makeblastdb -in $folder/$folder.fasta -dbtype nucl -out ../blast_db/$folder/$folder; done`

Server
`ls | while read folder; do /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/makeblastdb -in $folder/$folder.fasta -dbtype nucl -out ../blast_db/$folder/$folder; done`

`ls | while read folder; do mkdir ../05_blast_db_combined/$folder; ls | grep -v $folder | while read microbe; do cat $microbe/$microbe.fna >> ../05_blast_db_combined/$microbe/$microbe.fasta; done; /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/makeblastdb -in $folder/$folder.fasta -dbtype nucl -out ../blast_db/$folder/$folder; done`


.
ls 03_copy_of_02 | while read folder; do mkdir 05_blast_db_combined/$folder; ls done

ls | while read folder; do ls | grep


One thing to notice that makedb here i made it before realising there were problems with the probes so the db of those problematic microbes that cant produce the probes are still here, which is fine I guess???

7. Move the probes to another folder\
  `ls | while read folder; do mkdir ../../probes/$folder; done`

Filter the file, and wrote it as a fasta file cause blastn want it fasta
  `ls | while read folder; do cat $folder/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | cut -f4 | awk '{print ">\n"$0}' > ../probes/$folder/$folder.fasta; done`

  ls | while read folder; do mkdir ../04_probes/$folder; cat $folder/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | cut -f4 | awk '{print ">\n"$0}' > ../04_probes/$folder/$folder.fasta; done

8. Run blastn

`ls | while read folder; do mkdir ../../blastn/$folder; done`

Run blastn on the db, probes directory 

On the server it took about 4 hour for 58 samples, take it into account when writing time
`ls | while read folder; do blastn -db ../blast_db/$folder/$folder -query $folder/$folder.fasta -word_size 15 -ungapped > ../blastn/$folder/$folder.txt`



9. Filter uniq probes
Run quite slow, run in blastn folder
`ls | while read folder; do grep "Sbjct" -B 2 $folder/$folder.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - ../probes/$folder/$folder.fasta | grep -v ">" > ../uniq_probes/uniq_$folder.txt; done`
   
10. Try to find the mRNA location or whatever it is?

input requires three files: the original gtf, the paintshop result, and the uniq probe generate from 9.

10.1 Find the location of the probe

`grep -f probe_ap.txt paintshop_ap.tsv  | cut -f 1,2,3 > ../02_intermediate/location.txt`


`grep -f probes.txt paintshop_fc.tsv  | cut -f 1,2,3 > ../../02_intermediate/fc/location.txt`

10.2 Parse the important info from gtf file
`awk -F'\t' '$3 == "CDS"' ap.gtf | cut -f1,4,5,9 | h` 


FOR LOOP 10.1

`cat finish_microbe.txt | while read name; do mkdir find_gene_id_2/$name; grep -f uniq_probes/uniq_$name.txt extract_data/$name/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | cut -f 1,2,3 | sort -t$'\t' -k1,1 -k2,2n > find_gene_id_sorted/$name/probe_location_$name.txt`

alcanivorax_profundi
campylobacter_troglodytis
methylobacterium_isbiliense
streptococcus_phocae

`cat finish_unfinish_microbe.txt | while read name; do mkdir grep -f uniq_probes/uniq_$name.txt unfishish_extract_data/$name_*Balance.tsv | cut -f 1,2,3 | cut -d ":" -f 2 | sort -t$'\t' -k1,1 -k2,2n > find_gene_id_sorted/$name/probe_location_$name.txt; done`


FOR LOOP 10.2, run in find_gene_id_sorted folder

`ls | while read name; do awk -F'\t' '$3 == "CDS"' ../extract_data/$name/$name.gtf | cut -f1,2,3,4,5,6,7,8,9,10,11 | sort -t$'\t' -k1,1 -k2,2n  > $name/filter_$name.gtf; done`


`ls | while read name; do python3 ../../02_scripts/copilot_extract_gene_id.py $name/filter_$name.gtf $name/probe_location_$name.txt $name/gene_id_$name.txt; done`

### Do that again

`cat finish_microbe.txt | while read name; do mkdir find_gene_id_GENEID/$name; grep -f uniq_probes/uniq_$name.txt extract_data/$name/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | cut -f 1,2,3 > find_gene_id_GENEID/$name/probe_location_$name.txt`


`cat finish_unfinish_microbe.txt | while read name; do mkdir find_gene_id_GENEID/$name; grep -f uniq_probes/uniq_$name.txt unfishish_extract_data/$name_*Balance.tsv | cut -f 1,2,3 > find_gene_id_GENEID/$name/probe_location_$name.txt; done`
###

Find the universal probes
- Through blastn result, find all the matches and remove the none matches (opposite to what we did before)
`ls | while read folder; do grep "Sbjct" -B 2 $folder/$folder.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -f - ../probes/$folder/$folder.fasta | grep -v ">" > ../univ_probes/univ_$folder.txt; done`
- 


Bracken report file

ls | while read folder; do bracken -d /sw/data/Kraken2_data/prebuilt/k2_pluspf_20221209/ -i $folder/$folder.report -l S -r 100 -t 1 -o ../07_bracken_2/1_threshold/$folder.bracken_1; done

ls | while read folder; do bracken -d /sw/data/Kraken2_data/prebuilt/k2_pluspf_20221209/ -i $folder/$folder.report -l S -r 100 -t 8 -o ../07_bracken_2/8_threshold/$folder.bracken_8; done


cut -f1 * | sort | uniq > 2_uniq_speices.txt


# Univ probes
## makedb individual species
`cat finish_microbe.txt | while read microbe; do mkdir blast_db_individual/$microbe; /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/makeblastdb -in extract_data/$microbe/$microbe.fna -dbtype nucl -out blast_db_individual/$microbe/$microbe; done`

## cd 