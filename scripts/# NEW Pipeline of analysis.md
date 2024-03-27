# NEW Pipeline of analysis

# 1. Classification of microbes
[Metagenomics data](https://biorg.cs.fiu.edu/Smoking/) (RawData-CamposEtAl.zip) was taken from [this paper](https://www.microbiologyresearch.org/content/journal/acmi/10.1099/acmi.0.000497.v3#R52). The classification of microbes is done using the same pipeline from [my previous project](https://github.com/npxhuy/microbiome).

The classification step was done on UPPMAX's server.

With Bracken (v2.8) analysis, the result from Kraken were specified that reads be quantified at the taxonomic level of species level and only speices supported by 10 or more reads be counted.

Directory of working on UPPMAX: /proj/naiss2023-22-412/projects/microbiome/working/Hy/lungs.

Summary:
- Input: metagenomic data
- Output: a list of classified microbes
- Tool: FastQC, Trimmomatic, Kraken2, Bracken
- Script: check my github
- Other note: performed on UPPMAX. Have to filter the Bracken result (grep the first column of every result, uniq them, remove Homo sapiens as it was the only metazoan reference when running Kraken2)

# 2. Download data of microbes.
The number of classified speices with the threshold of 10 are 153. Using [Datasets command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) to download the FASTA and GTF file of the microbes.

The downloaded data is annotated, has a complete assembly level, comes from the RefSeq source, and is the latest assembly version. Some species have not been unclassified at species level yet so we exclude them. In the end we have 97 species left.

Summary:
- Input: a list of classified microbes
- Output: data of classified microbes
- Tool: Datasets command-line tools
- Script: *download_data.sh*
- Other note: perform locally and upload data to *Lunarc's server*

# 3. Probes design
## 3.1 Clean & Soft-link data
GTF files has # in the first few lines, in order for PaintSHOP to run, we need to clean them

Example code:
`ls | while read folder; cd $folder/ncbi_dataset/data/*/; do grep -v "#" genomic.gtf > genomic.tmp && mv genomic.tmp genomic.gtf; cd ../../../..; done`

The folder and file name of the reference genome data are named as their ID, so it's not a very good way to distinguish between files. Thus softlink them (*ln -s*) to another folder with different naming system would be better.

Example code:
`ls| while read microbes; do cd $microbes; ln -s ../../01_unzip_data/$microbes/ncbi_dataset/data/*/*.fna $microbes.fna; ln -s ../../01_unzip_data/$microbes/ncbi_dataset/data/*/*.gtf $microbes.gtf; cd ..; done`

Summary:
- Input: data of classified microbes
- Output: clean data of classified microbes
- Other note: perform on *Lunarc's server*

## 3.2 PaintSHOP
[PaintSHOP_pipeline](https://github.com/beliveau-lab/PaintSHOP_pipeline) was used to perform probe design. Follow the instruction on their github to run the pipeline, the pipeline was run on a loop repeatedly on all microbes' data.

The melting temperature, length of probes, and other settings need to be stated in the config.yml file. See the scripts for more info. Script stated below are an example that run only for my current working directory.

Summary:
- Input: clean data of classified microbes (fna + gtf)
- Output: intial probes
- Tool: PaintSHOP_pipeline
- Script: *paintshop_loop.sh*
- Other note: performed on *Lunarc's server*

# 4. Filtering
The large number of probes design by PaintSHOP could potentially off-binding (Probes designed for sp. A could potentially bind to sp. B). To advoid it [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) was used to remove all the potential off-binding probes.

Summary:
- Input: probes designed by PaintSHOP
- Output: "unique" probes
- Tool: BLAST+ (makeblastdb + blastn)
- Script: see example code
- Other note: performed on *Lunarc's server*, steps were stated below

## 4.1 Extract the probes' results and change them to fasta format
Because the PaintSHOP's results has so many different info/columns and we want the probes' column only for this part. Change to fasta format for BLAST+ to run normally.

Example code: `ls | while read folder; do mkdir ../04_probes/$folder; cat $folder/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | cut -f4 | awk '{print ">\n"$0}' > ../04_probes/$folder/$folder.fasta; done`



## 4.2 BLAST+: makeblastdb
Make database for running BLAST+. We will run blastn, for example, we have 3 sp. A,B,C; we will run sp. A's probes against B+C database. So the database for sp. A will include every other sp. exclude sp. A itself.

### 4.2.1 Combine fna file
Combine all the fasta file from every species excpet from the specie that will be run blastn against on.

Example code:
`ls | while read folder; do mkdir ..05_blast_plus/01_fasta_combined/$folder; ls | grep -v $folder | while read microbe; do cat $microbe/$microbe.fna >> ../05_blast_plus/01_fasta_combined/$folder/$folder.fasta; done; done`

### 4.2.2 makeblastdb
Run makeblastdb. -in take the fasta file (that we just combined on 4.2.1) -dbtype is database type, -out indicate the directory of output.

Example code:
`ls | while read folder; do mkdir ../02_makedb_01/$folder; /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/makeblastdb -in $folder/$folder.fasta -dbtype nucl -out ../02_makedb_01/$folder/$folder; done`

## 4.3: BLAST+: blastn
Run blastn. -db is the location of the database we made on 4.2.2, -query is the probe we made on 4.1.

Example code:
`ls | while read folder; do /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/blastn -db ../05_blast_plus/02_makedb_01/$folder/$folder -query $folder/$folder.fasta -word_size 15 -ungapped > ../05_blast_plus/03_blastn/$folder/$folder.txt; done`

## 4.4: Obtain "unique" probes:
Grep -v the result of blastn in 4.3 vs the original probe to have the unique probe

Example code:
`ls | while read folder; do grep "Sbjct" -B 2 $folder/$folder.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - ../../04_probes/01_original/$folder/$folder.fasta | grep -v ">" > ../../04_probes/02_unique/$folder/$folder.txt; done`

# 5. gene_id (& GeneID) finding of the probe

Input requires two files: 
1. the filtered gtf
2. the probe's location

## 5.1 Filtered GTF
On 97 samples of microbes, only 63 have their *db_xref* (linking DNA sequence records to other external databases) to Entrez Gene Database. We aim to select probes that will land closely to each other (on similar region or having probe having same gene_id)

To seperated between two kind of gtf file, two lists were made including gene_id_list.txt (gtf files do not have GeneID section) and GeneID_list.txt (gtf files that have GeneID section)

Filter it to have the chromosome name, gene location (start and stop position) and gene id.

Example code:
- For GeneID files:

`grep "GeneID" -m 1 03_copy_of_02/*/*.gtf | cut -d "/" -f2 > GeneID_list.txt`

`cat GeneID_list.txt | while read name; do mkdir 06_gene_id/$name; grep "GeneID" 03_copy_of_02/$name/$name.gtf | cut -f1,4,5,9 | sort -t$'\t' -k1,1 -k2,2n > 06_gene_id/$name/${name}_filter.gtf; done`
- For gene_id files:

`grep "GeneID" -m 1 03_copy_of_02/*/*.gtf | cut -d "/" -f2| grep -v -f - microbes_list.txt > gene_id_list.txt`

`cat gene_id_list.txt| while read name; do mkdir 06_gene_id/$name; awk -F'\t' '$3 == "CDS"' 03_copy_of_02/$name/$name.gtf | cut -f1,4,5,9 | sort -t$'\t' -k1,1 -k2,2n  > 06_gene_id/$name/${name}_filter.gtf; done`


awk -F'\t' '$3 == "CDS" || $3 == "transcript"' 03_copy_of_02/$name/$name.gtf | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F';' '{print $1"; "$2";"}' | sort -t$'\t' -k1,1 -k2,2n  > 06_gene_id/$name/${name}_filter.gtf; done

THIS CODE
awk -F'\t' '$3 == "CDS" || $3 == "transcript"' ../03_copy_of_02/Acinetobacter_baumannii/Acinetobacter_baumannii.gtf | cut -f 1,4,5,9 | cut -d ";" -f 1,2 | sed 's/gene_id\|transcript_id\|;\|\"//g'


## 5.2 The probes' location
Obtaining by parse the paintshop result and the uniq probe generate from **4.4**

Example code:

`cat microbes_list.txt | while read name; do grep -f 04_probes/02_unique/$name/$name.txt 03_copy_of_02/$name/pipeline_output/03_output_files/01_dna_probes/*Balance.tsv | sort -t$'\t' -k1,1 -k2,2n > 06_gene_id/$name/${name}_probe_location.txt; done`

## 5.3 Find the gene_id from those two files

`ls | while read name; do python3 ../../02_scripts/copilot_extract_gene_id.py $name/${name}_filter.gtf $name/${name}_probe_location.txt $name/${name}_gene_id.txt; done`


# 6. Readout
## 6.1 Readout sequences
### 6.1.1 blastn the readout sequence
From [MERFISH's paper 2016](https://www.pnas.org/doi/full/10.1073/pnas.1612826113#sec-2) SI (pnas.1612826113.sd01.xlsx). I found that they use 16 readout sequences in their probes. "The 20-nt, three-letter readout sequences were designed by generating a random set of sequences with the per-base probability of 25% for A, 25% for T, and 50% for G."

Make individual database for BLAST, so we can run blastn of the 16 readout sequences against the database to see if they have any homology. 

"subset of these sequences with no cross-homology regions longer than 11 contiguous bases." => So the wordsize will be 11.

Example code:
`ls | while read folder; do /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/makeblastdb -in ../../03_copy_of_02/$folder/$folder.fna -dbtype nucl -out $folder/$folder; done`

Turned out they are most of the query matched with the database hmm. Still we need to extract the ID and the location of the matches to evaluate.

### 6.1.2 Preparation for evaluation of the readout sequence's location
Making the input files, including the probes' location file and readouts' location file.

Example code:
`cat microbes_list.txt | while read microbe; do /home/npxhuy/04_tools/ncbi-blast-2.15.0+/bin/blastn -db 05_blast_plus/04_makedb_individual/$microbe/$microbe -query 07_readout/01_seq/readout_sequences.txt -word_size 11 -ungapped > 05_blast_plus/05_blastn_readout_sequence/$microbe/$microbe.txt; done`

We then need the probes' location file to compare, we have it from 5.2, but we gonna extend the probes' location cause this readout sequence will be attach on either side of the probes, to be precise we want to know the 20 nucleotides before and after the probe loction if it matchs with the readout sequence or not, but to simplify and to boarden the range, we extend it to 60.

Example code: `awk '{$2=$2-60; $3=$3+60; print}' filename`

### 6.1.3 Evaluation of the readout sequence's location

Now we have two files needed, we want to run the matching scripts to see if the readout sequence may match with the upper and lower region of the probe or not.

Scripts: probe_readout_matching.py

Example run:
`ls | while read folder; do python3 /home/npxhuy/02_scripts/probe_readout_matching.py $folder/${folder}_probe_location_extention.txt $folder/${folder}_readout_seq_location.txt $folder/${folder}_matching.txt`

cut -f 1,2,3,4,5,6,7,8,9,10 Acinetobacter_baumannii_probe_location.txt > b.txt
python3 gene.py Acinetobacter_baumannii_filter.gtf b.txt a.txt






/home/npxhuy/lu2023-17-27/hy/thesis/03_data/03_copy_of_02/Acinetobacter_baumannii/pipeline_output/03_output_files


cat /home/npxhuy/lu2023-17-27/hy/thesis/03_data/03_copy_of_02/Acinetobacter_baumannii/pipeline_output/03_output_files/01_dna_probes/Acinetobacter_baumannii_all_newBalance.tsv

/home/npxhuy/lu2023-17-27/hy/thesis/03_data

TTGATCATGGCTCAGATTGAACGCTGGCGG

/home/npxhuy/lu2023-17-27/hy/thesis/04_tools/ncbi-blast-2.15.0+/bin/blastn -db /home/npxhuy/lu2023-17-27/hy/thesis//03_data/05_blast_plus/02_makedb_01/Acinetobacter_baumannii/Acinetobacter_baumannii -query new_rna.txt -word_size 15 -ungapped > blastn_result.txt


grep "Sbjct" -B 2 blastn_result.txt | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - new_rna.txt | grep -v ">" > final_rna_probe.txt