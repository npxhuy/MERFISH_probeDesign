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

