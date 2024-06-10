# Pipeline of analysis
Note: some tasks are more computationally intense than others and more complicated; their location is in the *scripts* folder. If 
tasks/scripts are lightweight example code snippet can be found in this document except were noted. Change the example code to fit your file's location/directory/naming.

The figure below show the pipeline of analysis with the rectangle boxes equivalent to the sections in this README.
![alt text](https://github.com/npxhuy/MERFISH_probeDesign/blob/main/figure/flowchart_MD.png)

## 1. Classification of Microbes and Abundance estimation.
[Metagenomics data](https://biorg.cs.fiu.edu/Smoking/) (RawData-CamposEtAl.zip) was taken from [this paper](https://www.microbiologyresearch.org/content/journal/acmi/10.1099/acmi.0.000497.v3#R52). The classification of microbes is done using the same pipeline from [my previous project](https://github.com/npxhuy/microbiome). Please refer to that GitHub link to do this step.


For short
1. You must QC the metagenomic data with FastQC for decision-making. You can skip the trimming step if the data gives good results with FastQC.


2. You will need to run Kraken2 on the trimmed data. You can always make your database if you do not want to use the PlusPF database. To build a database, please refer to Kraken2 [manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown).


3. You will need to run Bracken on the Kraken2 results. If you make your database for Bracken, you need to build the Bracken k-mer database within your database. Please refer to Bracken [manual](https://github.com/jenniferlu717/Bracken)\
Having the classification of microbes, the first column of every bracken report is the name of species, cut them out and put them in a file, which will be used later on.


## 2. Raw data processing
Now that you have the list of microbes, you can download their genome from NCBI, for example. If you can find an alternative way to download microbe data, it's good for you; if you already have them, it's also good for you. If not, use [datasets command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) to download the FASTA and GTF file of the microbes from NCBI. For this research, we used PaintSHOP to design target probes using FASTA and GTF files. Example code for this task is *download_data.sh*. You will need the microbes list from **1** for this task.


The genome files are in zip format; unzip them! GTF files have bunches of lines starting with #, which can prevent PaintSHOP from running; remove it.


Example code:
`ls | while read folder; cd $folder/ncbi_dataset/data/*/; do grep -v "#" genomic.gtf > genomic.tmp && mv genomic.tmp genomic.gtf; cd ../../../..; done`


The folder and naming are pretty tricky because of multiple folders and using accession numbers instead of scientific names for the FASTA and GTF files, which were all named genomic.gtf. To make it easier, you can soft-link them to another folder.


Example code:
`ls| while read microbes; do cd $microbes; ln -s ../../01_unzip_data/$microbes/ncbi_dataset/data/*/*.fna $microbes.fna; ln -s ../../01_unzip_data/$microbes/ncbi_dataset/data/*/*.gtf $microbes.gtf; cd ..; done`


## 3. Probe design
### 3.1 Target sequence
#### 3.1.1 Probe design
[PaintSHOP_pipeline](https://github.com/beliveau-lab/PaintSHOP_pipeline) was used to perform probe design. Follow the instructions on their GitHub to run the pipeline. The pipeline was run on a loop repeatedly on all microbes' data.


The melting temperature, length of probes, and other settings need to be stated in the config.yml file. See the script *paintshop_loop.sh* for more info. This is a modified script base on their *run_pipeline.sh* script, to make it into a loop and write config.yml file for every species, that's the fundamental of the scripts, in case it's not working with you device, I suffered it too.

 The error may occur when you give a high Tm with a short probe length, making the tool unable to produce any probes. With that Tm and length, some probes can satisfy the input. You can consider separating the code into two parts: generating a config.yml file for all the microbes and running PaintSHOP on the data for a less complicated task.


Example of output:
> NZ_CP043953.1   174     203     ATTGGCTGAACAATTGTCAGAAGGGCGGGT  42.480  100.000 0.000   0       0.307   0       +\
NZ_CP043953.1   256     285     AGTGAACAGCCTGCAACAACTACAGCAGCT  42.260  96.458  0.000   0       0.285   0       +\
NZ_CP043953.1   417     446     TGAAGGCCGTTCTAACCAAATGGCAGCAGA  42.200  100.000 0.000   0       0.265   0       +




#### 3.1.2 Make databases for blastn


The probe was screened against the human genome and with other microbes with [Blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). To run Blast+, we need a blast database of that genome.


Download a version that fits your current platform, and we can start with the making of the database using makeblastdb tool. For this step, two main db will be made: 1. human genome, 2. microbes genome.


1. Human genome database\
Download the fasta file from [here](https://www.ncbi.nlm.nih.gov/genome/guide/human/). Extract the file if needed and run *makeblastdb* on it.\
Example code: `makeblastdb -in $human.fasta -dbtype nucl -out $human_db`


2. Microbe genome database\
To run screen one speice to all other species other than that species, we need to combine the data of other species and then *makeblastdb*
- Example code of combine fasta file: check *combine_fasta_for_makeblastdb.sh*
- Example code of *makeblastdb* on the combined fasta file: check *makeblastdb_microbe.sh*


#### 3.1.3 Screening with blastn + filtering
This step will be repeated many times. When we have the database, we run blastn to screen the probe on the database. Having the blastn result, we parse it, extract the found sequences (which means those sequences are not unique) and filter it back to the original file.


1. Make the probes result in fasta format.\
Example code: `cat $probe_result.tsv | cut -f4 | awk '{print ">\n"$0}' > $probe_result.fasta`\
Example output:
> \> \
ATTGGCTGAACAATTGTCAGAAGGGCGGGT\
\> \
AGTGAACAGCCTGCAACAACTACAGCAGCT\
 \> \
TGAAGGCCGTTCTAACCAAATGGCAGCAGA
2. Run blastn. This could be heavy in term of computational resources, but not at this step yet.\
Example code:`blastn -db $db_name -query $probe_result.fasta -word_size $number -ungapped > $blast_result_file`
3. Parse blast result.\
Example code: `grep "Sbjct" -B 2 $blast_result_file | grep "Query" | awk '{print $3}' | sort | uniq | grep -v -f - $probe_result.fasta | grep -v ">" > $unique_probe`
4. (Optional if you want to repeat the process) Change the $unique_probe to fasta format.\
Example code: `awk '{print ">\n"$0}' $unique_probe > $unique_probe.fasta`


Run these steps to screen the human genome, obtaining the probes that are unique to the human genome. Screen the unique probes on the human genome on microbes genome, obtaining the unique probes on human + microbes.


In step 2, technically, it can run normally with a single thread, taking around 5-10 mins without using any additional core power. However, you can utilise multiprocessing for faster time processing, which will be helpful in the Readout Sequence design step. Example code for multi-threading will be provided for that step later on.


#### 3.1.4 Filter to RNA probe.
This step required two inputs: 1. a filtered gtf file, and 2. a probes' location file
1. Make filtered gtf files.\
Example code: `awk -F'\t' '$3 == "CDS" || $3 == "transcript"' $microbe.gtf | cut -f 1,4,5,9 | cut -d ";" -f 1,2 | sed 's/gene_id\|transcript_id\|;\|\"//g' > $microbe_filter.gtf`\
Example output:\
> NZ_CP043953.1   1       1395     F3P16_RS00005  unassigned_transcript_1\
NZ_CP043953.1   1496    2641     F3P16_RS00010  unassigned_transcript_2\
NZ_CP043953.1   2659    3738     F3P16_RS00015  unassigned_transcript_3
2. Make probe location files. Obtaining by parsing the PaintSHOP result and the unique probes generated from **3.1.1** and the unique probe from **3.1.3**\
Example code: `grep -f $unique_probe $probe_result.tsv | sort -t$'\t' -k1,1 -k2,2n | cut -f 1,2,3,4,5,6,7,8,10,11 > $probe_location.txt`
Example output:\
> NZ_CP043953.1   4164    4193    AAGTTTCAGGTGGCTTACACGGCGTAGGTG  42.210  100.000 0.000   0       0       +\
NZ_CP043953.1   7868    7897    TTTATCGACTTTAGCGGGAGTTGGAGCGGC  42.100  100.000 0.000   0       0       +\
NZ_CP043953.1   9963    9992    GCTCTAGTAGATCGTAGGCTGCGTGAGGGT  42.080  100.000 0.000   0       0       +
3. Parse from those two files to get the result. Script *extract_gene_id.py* was written for this step.
Example code: `python3 extract_gene_id.py $microbe_filter.gtf $probe_location.txt $output_gene_id.txt`\
Example output: No_Match indicates that the probe won't land in the coding region.
> Chr     Start   Stop    Seq     Tm      on-targ off-targ        repeat  max k-m strand  GeneID  Transcript\
NZ_CP043953.1   4164    4193    AAGTTTCAGGTGGCTTACACGGCGTAGGTG  42.21   100.0   0.0     0       0       +        F3P16_RS00020  unassigned_transcript_4 \
NZ_CP043953.1   20214   20243   AAGCTGTACGGTGCTTAAGTGCACAGTGCT  42.2    98.981  494.522 0       3       +       No_Match        No_Match


For each species, we take randomly a maximum of 50 probes and a minimum of 20; any species that doesn't have at least 20 probes was discarded. Probe that has "No_Match" and off-target score greater than 0.0 will be removed. Using the code *selecting_probe.sh*. This script was run in the directory where having all the subdirectorys' name are the microbes' name. In each subdirectory the *$output_gene_id.txt* from this step was named with the microbe name + _ + gene_id.txt. If you use a different naming or directory, you should chagne the *selecting_probe.sh* a little bit.

Maybe not every species may went through all the extreme filtering steps, let's see how many we have left.

Examle code: `wc -l */*gene_id.txt | awk '{print $2}' | cut -d / -f 1 | head -n -1 > final_microbes.txt`

This will be our final list of microbes!

### 3.2 Readout sequence design
Readout sequence design started with generating all possible combinations of A, T, and G with a % of ~25,25,50, respectively. The % varies a little bit; for example, for a 20-nt sequence, I could have 6A5T9G. This aims to generate as much sequence as possible for the readout selection. Use the script *generate.py* for this process.


Example usage: `python3 generate.py --output $sequence --target A=5,G=10,T=5 --constraint GGGG=0 --process 8 --verbose`


This code can be generated quite effectively using the multicore processing within Lunarc's server. With only a three-word combination, it's moderately fast, taking 15-20 mins with four threads and less than 5 mins with eight threads. However, when you want to generate sequences with all ATCG and a length of 20, I warn you it might take forever. This code runs well for a 3-base combination, but for 4, I do not think so.


Another code was used to filter down the generated sequences. Because they are usually one base different from each other, running BLAST+ on this sequence is not a good idea for one sequence with a homology of 11. There is a high chance that the following thousand sequences also share the homology. Therefore, another code, *filter.py*, was used to filter it down first. It randomly separates the sequence into different groups and compares sequences within the group to each other in multiple rounds.


Example usage: `python3 filter.py --input $sequence --output $filter_sequence --distance 4 --iteration 5 --chunk 1000 --process 8 --verbose`


After this step, repeat part **3.1.3** for the rest of the sequence (word size is 11 for human transcriptome and 12 for microbes), and screen the potential readout sequences with themselves so they are different enough from each other. Select 16 after filtering. Readout probes are the reverse complement of readout sequences.

**OPTIONAL**\
The part of repeating **3.1.3** could be very time-consuming, as there are millions of potential readouts that need to be filtered down with BLAST+. Hence, check the script *example_for_multiprocessing.py* to see how I did it. I recommend you to follow it. Or even improve it in your own way.



### 3.3 Primers design


I started with this file, a [FASTA file of the 240K 25mer barcodes](http://elledgelab.med.harvard.edu/wp-content/uploads/2013/07/bc25mer.240k.fasta_.zip) from Elledge lab, Harvard. 25kmers was trimmed to 6 20kmers sequences using *cutting25mer.py* (the first one start from 1-20, and then 2-21,...).\
Example usage: `python3 cutting25mer.py bc25mer.240k.fasta cut_20mer.fasta`


Those sequences then went through primer filtering (GC content, Tm, hairpin formation, etc.) using *primerFilter.py*.\
Example usage: `python3 primerFilter.py cut_20mer.fasta left_over.fasta`


After this, repeat part **3.1.3**. Word size is 11 for human transcriptome and 12 for microbes. After this select two sequences being the potential primers.


## 4. Codebook design
The codebook was adapted from a [study on MERFISH](https://www.science.org/doi/10.1126/science.aaa6090) (in supplementary part). They provided 16-bit MHD4 code (140 codes). You can design the code book using [La Jolla Covering Repository Tables](https://ljcr.dmgordon.org/cover/table.html). With hamming distance 4, hamming weight 4, and 16 bit, you should select t=3,k=4,v=16. Increasing the bit will increase the v.


With a table of C(16,4,3), we will have the following information:
> 1  2  3  4\
1  2  5  6\
1  2  7  8


These numbers indicate the position of 1 in a binary barcode, so for the first line, the barcode should be 1111000000000000. We first need to download the HTML file (website) with the position to generate the binary barcode. Right-click -> save as - so we can save the page's html file. Use *binary_barcode_generation.py* to generate the binary barcode.


Example usage: `python3 binary_barcode_generation.py $html_file $output`

## 5. Combine data
We now have
1. Target probes (from 3.1)
2. Readout sequences (from 3.2)
3. Primers (from 3.3)
4. Codebook (from 4)
5. Final microbes list (from 1,2 and the last step of 3.1)

We can start combining them in the result table.

### Making intermediate files from input
1. 
Create a directory called *01_codebook*, copy needed files stated below in.\
We first will combine **Final list** with **Codebook**

Microbes list
> Acinetobacter_baumannii\
Actinopolyspora_erythraea\
Agathobacter_rectalis\
Anaerococcus_mediterraneensis

Codebook
> 1111000000000000\
1100110000000000\
1100001100000000\
1100000011000000

Run the script *assign_barcode.py* to combine them into an intermediate file. It will randomly assigne a barcode to a species, any barcode left was left with NO_ASSIGNED
Example usage: `python3 assign_barcode.py codebook.txt microbes_list.txt assigned_output.txt`
> Acinetobacter_baumannii	[4 6 12 14]	0001010000010100\
Actinopolyspora_erythraea	[1 6 9 14]	1000010010000100\
Agathobacter_rectalis	[11 12 13 14]	0000000000111100\
Anaerococcus_mediterraneensis	[5 7 10 12]	0000101001010000

2. 
Create a directory called *02_readout_seq*, copy needed file stated below in.\
16 readout sequences with their own round number (from 1 to 16)
Example code: `nl -n ln -w 1 -s ' ' readout.txt > temp; rm readout.txt; mv temp readout.txt`

3. 
Create a directory called *03_selected_probes* and copy the result from **3.1** in, you might want to change the name of the file to without _gene_id.txt, at least in my case.
### Combining them into one table
Now that you got everything ready, a script *combine.py* was run outside the three above directories, and it will produce a file call *probe_table.tsv*
You have to open the code, edit line 76-77 with your desire primers. I didnt do it smarter because I'm sort on time with my thesis, not that I do not have coding skills.

Example output:
>Specie  Chr     Primer_1        Readout_1       Readout_2       Target  Readout_3       Readout_4       Primer_2\
Acinetobacter_baumannii NZ_CP043953.1   CGTGTTAGTGGCCCGGGTCT            AGATGGATAGGGTTATGAGT    AAAGAAGCCAGACAATTCGTCCCGGATGGC  AGATGTAGGTTAGGTGAGAG    GTTAGAGAGAGAGAGGGTTT    GGCCGCGACTAGGTAAGCCT\
Acinetobacter_baumannii NZ_CP043953.1   CGTGTTAGTGGCCCGGGTCT    AGATTAGGAGGGTAGTTATG            ACATCCGGAGTCGAGGTTTCGCTTACACCC  AGATGTAGGTTAGGTGAGAG    GTTAGAGAGAGAGAGGGTTT    GGCCGCGACTAGGTAAGCCT\
Acinetobacter_baumannii NZ_CP043953.1   CGTGTTAGTGGCCCGGGTCT    AGATTAGGAGGGTAGTTATG    AGATGGATAGGGTTATGAGT    TGTTCTTAGTCTCGTGTTAGGGTCCGGGCT          GTTAGAGAGAGAGAGGGTTT    GGCCGCGACTAGGTAAGCCT


### Well done following the pipeline. Thanks for reading!!!






















