List of rules
# References files

### 1. rule *parse_genome*
- input: FASTA file
- output: multiple output files, located in the 01_reference_files/01_chorm_names +
and a FASTA file in 01_reference_files/02_multi_fasta\
- params: mostly the same in other rules, include mfree (the memory for this run) and h_rt(hour of run time)
- scripts: parse_genome.py
*Summary*: produce references file, chromosome name file and multifasta file are two important files here.
1. **Import necessary module**: The script starts by importing the `SeqIO` module from the `Bio` package.
    ```python
    from Bio import SeqIO
    ```

2. **Configure file paths and setup variables**: The script sets up file paths and other configuration options using the `snakemake` workflow management system. It also sets up lists to keep track of sequences that are removed based on different criteria.
    ```python
    INPUT_FA = snakemake.input.fasta
    OUTPUT_FA = snakemake.output.multi_fasta
    CHROM_NAMES = snakemake.output.chrom_names
    ...
    EXCL_CHROMS = []
    if 'exclude_chroms' in snakemake.config:
        EXCL_CHROMS = snakemake.config['exclude_chroms']
    ```

3. **Define the main function**: This function does the bulk of the work.
    ```python
    def main():
    ```

4. **Load the FASTA file**: The script reads in the FASTA file using `SeqIO.parse`.
    ```python
    recs = SeqIO.parse(snakemake.input[0], 'fasta')
    ```

5. **Loop over each sequence record**: For each sequence, it checks the sequence identifier (`seq.id`) to determine whether the sequence should be included in the final output or filtered out.
    ```python
    for seq in recs:
    ```

6. **Filter sequences based on their identifiers**: Sequences are filtered out based on whether their identifier contains certain substrings ('Un_', '_random', '_alt', '_hap', '_fix') or is in the `EXCL_CHROMS` list. If a sequence is filtered out, its identifier is added to the appropriate list (e.g., `removed_unk`, `removed_rand`, etc.).
    ```python
    if 'Un_' in seq.id:
        removed_unk.append(seq.id)  # keep a record of this chromosome
        filtered.append(seq)        # add to multifasta for bowtie2 index
    ...
    ```

7. **Write the filtered sequences to a new FASTA file**: After all sequences have been processed, the script writes the `filtered` sequences to a new FASTA file.
    ```python
    with open(OUTPUT_FA, 'w') as outfile:
        SeqIO.write(filtered, outfile, 'fasta')
    ```

8. **Write the identifiers of the filtered sequences and the identifiers of the included sequences to separate output files**: It also writes the identifiers of the filtered sequences and the identifiers of the included sequences to separate output files.
    ```python
    with open(CHROM_NAMES, 'w') as outfile:
        outfile.write('\n'.join(chrom_names))
    ...
    ```

9. **Print a success message and the total runtime**: Finally, the script prints a success message and the total runtime.
    ```python
    print('DONE!')
    print('Script took {} seconds'.format(time.time() - start))
    ```

### 2. checkpoint *split_chroms*
- input: multi_fasta and chrom_names files from **1.**
- output: chrom_fasta_dir
- scripts: split_chroms.py
*Summary*: split the multi fasta file into seperated files, each file contains the sequence of its chromosome.
1. **Import necessary modules and configure file paths**: The script starts by importing the `os` and `SeqIO` modules. It then sets up file paths using the `snakemake` workflow management system.
    ```python
    import os
    from Bio import SeqIO
    INPUT_FA = snakemake.input.fasta
    CHROM_NAMES = snakemake.input.chrom_names
    CHROM_FASTA_DIR = snakemake.output.chrom_fasta_dir
    ```

2. **Create output directory**: The script creates an output directory where the individual chromosome files will be saved.
    ```python
    os.makedirs(CHROM_FASTA_DIR, exist_ok=True)
    ```

3. **Define the main function**: This function does the bulk of the work.
    ```python
    def main():
    ```

4. **Load the genomic FASTA file**: The script reads in the FASTA file using `SeqIO.parse`.
    ```python
    recs = SeqIO.parse(INPUT_FA, 'fasta')
    ```

5. **Load the list of chromosome names**: The script calls the `get_chrom_names` function to get a list of chromosome names.
    ```python
    chrom_names = get_chrom_names()
    ```

6. **Save individual chromosome files**: The script loops over each sequence record in the FASTA file. For each sequence, it checks if the sequence identifier (`seq.id`) is in the list of chromosome names. If it is, the sequence is written to a new file in the output directory. If it's not, the sequence is skipped.
    ```python
    for seq in recs:
        if seq.id not in chrom_names:
            print('Skipping non-canonical chromosome: ', seq.id)
            continue
        file_path = os.path.join(CHROM_FASTA_DIR, f'{seq.id}.fa')
        with open(file_path, 'w') as outfile:
            print(f'Writing {seq.id}')
            SeqIO.write(seq, outfile, 'fasta')
    ```

7. **Define the `get_chrom_names` function**: This function reads in a file of chromosome names and returns a list of these names.
    ```python
    def get_chrom_names():
        with open(CHROM_NAMES, 'r') as infile:
            text = infile.read().strip(' \n')
        chrom_names = text.split('\n')
        return chrom_names
    ```

8. **Run the main function**: If the script is run directly (not imported as a module), the `main` function is called and the runtime is printed.
    ```python
    if __name__ == '__main__':
        import time
        start = time.time()
        main()
        print('Script took {} seconds'.format(time.time() - start))
    ```

### 3. checkpoint *parse_gtf*
- input: annotation gtf file and chrom_names from **1**
- output: a **bed** file and a directory of chromo_df_dir where it will have **dat** files for each chromosome make from **pickle** module
- scripts: parse_gtf.py\
*Summary*: parse the *gtf* file, filter it for exons, save the filtered data in a *BED* file. This BED file contains the exon annotations for each transcript in the input GTF file, with each row representing an exon and each column representing a different attribute of the exon (such as its start and end positions, the transcript it belongs to, its score, its strand, its transcript version, and its gene ID). For each chromosome there will be a *dat* file, not sure what they will do with that file
1. **Import necessary modules and configure file paths**: The script starts by importing the `os`, `pickle`, `numpy`, and `pandas` modules. It then sets up file paths using the `snakemake` workflow management system.
    ```python
    import os
    import pickle
    import numpy as np
    import pandas as pd
    CHROM_NAMES = snakemake.input.chrom_names
    GTF_PATH = snakemake.input.gtf
    BED_PATH = snakemake.output.bed
    CHROM_DF_DIR = snakemake.output.chrom_df_dir
    ```

2. **Create output directory**: The script creates an output directory where the individual chromosome files will be saved.
    ```python
    os.makedirs(CHROM_DF_DIR, exist_ok=True)
    ```

3. **Define the main function**: This function does the bulk of the work.
    ```python
    def main():
    ```

4. **Load the GTF file**: The script reads in the GTF file using `pandas.read_csv` and renames the columns for clarity.
    ```python
    df = load_gtf(GTF_PATH)
    ```

5. **Filter the GTF data**: The script filters the GTF data to only include exons on canonical chromosomes.
    ```python
    df = parse_and_filter(df)
    ```

6. **Split the GTF data into individual chromosome files**: The script loops over each unique chromosome in the GTF data. For each chromosome, it filters the GTF data to only include that chromosome and saves the filtered data to a file.
    ```python
    chroms = df['seqid'].unique()
    for chrom in chroms:
        chrom_df = df[df['seqid'] == chrom]
        df_path = os.path.join(CHROM_DF_DIR, f'{chrom}_filtered_gtf.dat')
        pickle.dump(chrom_df, open(df_path, 'wb'))
    ```

7. **Save the BED file**: The script saves a subset of the GTF data as a BED file.
    ```python
    bed_df = df[['seqid', 'start', 'end', 'transcript_id', 'score', 'strand', 'transcript_version', 'gene_id']]
    bed_df.to_csv(BED_PATH, sep='\t', index=False, header=None)
    ```

8. **Define the `load_gtf`, `parse_and_filter`, `get_chrom_names`, and `parse_attributes` functions**: These functions are used to load the GTF file, filter the GTF data, get the list of canonical chromosome names, and parse the attributes column of the GTF data, respectively.
    ```python
    def load_gtf(file_path):
    def parse_and_filter(df):
    def get_chrom_names():
    def parse_attributes(data):
    ```

9. **Run the main function**: If the script is run directly (not imported as a module), the `main` function is called.
    ```python
    if __name__ == '__main__':
        main()
    ```

### 4. rule *flatten_isoforms*
- input: the *dat* file and the filtered gtf file in **3**, and basically all files in **3** as *upstream* (but didnt see any in use of that input in the py script)
- output: an iso_flat **bed** file
- scripts: iso_flatten.py
*Summary*: Open the pickled chromosome dataframe (binary file), get flattened exons and produce the flatten isoform data. Flattening exons may involve combining or merging adjacent exons that belong to the same transcript or gene, resulting in a more concise representation. Still a bit unclear.\
The Python code consists of two functions: `flatten_isoforms` and `get_flattened_exons`.

1. `flatten_isoforms(df)`: This function takes a DataFrame `df` as input and performs the following steps:
    - Groups exons by gene: The DataFrame is grouped by the 'gene_id' column. For each group, the first 'seqid', the tuple of 'start' values, the tuple of 'end' values, the first 'score', the tuple of 'strand' values, and the tuple of 'transcript_version' values are aggregated.
    - Performs isoform flattening: The `get_flattened_exons` function is applied to each row of the grouped DataFrame. The result is a Series of lists of exon data. This Series is flattened into a single list of exon data, which is then converted into a DataFrame.
    - Returns the DataFrame: The DataFrame of flattened exon data is returned.

2. `get_flattened_exons(row)`: This function takes a row of a DataFrame as input and performs the following steps:
    - Converts input coordinates to numpy arrays: The start and end coordinates of the exons are converted to numpy arrays.
    - Calculates the total span of all exons: The total span of all exons is calculated by finding the minimum start coordinate and the maximum end coordinate.
    - Zero shifts coordinates: The start and end coordinates are shifted so that the minimum start coordinate is now 0.
    - Calculates exon coverage: The coverage of each exon, which is the number of isoforms that include that exon, is calculated.
    - Collapses isoforms as needed: If the maximum coverage is greater than 1, the function finds the start and end coordinates of the maximally shared segments and updates the output coordinates accordingly.
    - Formats output annotation data for shared segments: The output annotation data for the shared segments is formatted and appended to a list.
    - Returns the list of segments: The list of segments is returned.

### 5. def *get_flat_beds*
- *Summary*: Return a list of file paths from **3** and **4** that will be use for **6**

1. Constructs a path to the output files of the `parse_gtf` checkpoint. The `**wildcards` syntax is used to pass the values of all the wildcards as keyword arguments to the `get` method of the checkpoint. This allows the function to dynamically determine the output of the checkpoint based on the wildcards. The output directory of the checkpoint is joined with a filename pattern that includes a wildcard `{chrom_gtf}_filtered_gtf.dat`.
2. Uses the `expand` function to construct a list of paths to the output files of the `flatten_isoforms` rule. The `glob_wildcards` function is used to extract the values of the `chrom_gtf` wildcard from the wildcard path.
3. Returns the list of file paths. This list can be used as the input to a rule, allowing the rule to operate on all the files produced by the `flatten_isoforms` rule for each chromosome.

### 6. rule *merge_flat_annotations*
- input: the file path from **5**
- output: an iso flat bed file

### 7. rule *build_bowtie*
- input: the multi_fasta from **1**
- output: build a bowtie2 index from the filtered multi-fasta genome
- extra info: The command bowtie2-build is used to build this index from a reference genome. The resulting index files are used by Bowtie2 for the alignment of short read sequences. This is a necessary step before running Bowtie2 for sequence alignment.

### 8. rule *build_jellyfish*
- input: the multi_fasta from **1**
- output: build a jellyfish index from the filtered multi-fasta genome.
- *Summary*: counting the frequency of each unique sequence of 18 bases in the input, discarding any that appear less than 2 times, and writing the counts to the output file.
The `jellyfish count` command in the Snakemake rule is used to count k-mers in DNA sequences. Here's a breakdown of what this command does:

- `-m 18`: This option sets the k-mer size to 18. Jellyfish will count all unique sequences of 18 bases in the input.

- `-s 3300M`: This option sets the initial hash size to 3300M items. This is used to optimize the memory usage of Jellyfish.

- `-o {output}`: This option specifies the output file for the k-mer counts. The `{output}` is a placeholder that Snakemake will replace with the actual output file path.

- `--out-counter-len 1`: This option sets the counter length for the output. This is the number of bytes used to store the count of each k-mer.

- `-L 2`: This option tells Jellyfish to discard k-mers that appear less than 2 times in the input.

- `{input}`: This is the input file(s) for Jellyfish to count k-mers in. The `{input}` is a placeholder that Snakemake will replace with the actual input file path(s).\

# Intermediate files

### 9. rule *design_probes*
- input: the multi_fasta from **1**
- output: a fastq file contains inital probes for each chromosome
- script: blockParse_unmasked.py\ the script is so long i decided to look at it later, basically it will result in the initial probe

Possible explanation of the script is: screens for prohibited sequences such as homopolymeric runs and “N” bases and allows users to specify allowable ranges of probe length, percent G + C content (GC%), and adjusted Tm calculated by using nearest-neighbor thermodynamics

### 10. rule *format_exising_probes*
- optinal rule if we provide probe it will skip **9**

### 11. rule *align_probes*
- input: result from **9**
- output: a `bam` binary file
- *Summary*: align probes to reference genome using bowtie2 and samtools, but not sure what bowtie2 and samtools would do, will figure out later

### 12. rule *sam_to_pairwise*
- input: result from **11**
- output: a txt .out file
- *Summary*: extract duplex information from alignment results

### 13. rule *extract_alignment_scores*
- input: result from **11**
- output: a txt file showing alignment score
- *Summary*: extract alignment scores from pairwise alignment

### 14. rule *parse_pairwise*
- input: result from **12 and 13**
- output: a csv file showing alignment score
- Script: parse_pairwise.\
*Summary*: construct alignment dataframe using pairwise alignment info

### 15. rule *predict_duplex*
- input: result from **14** and a binary pickled_models file
- output: a csv file contain the duplexing
- Script: XGB_predict.py\
 *Summary*: generate predictions of the duplexing probability of probe candidates and their Bowtie2 alignments. predict duplexing probabilities using pre-trained XGBoost model

### 16. rule *score*
- input: result from **15**
- output: a bed file having the on target score and off target score
- Script: output_bed.py\
Not understand how they calcuate yet.

### 17. rule *filter_n_probes*
- input: result from **16**
- output: a bed file, filter from input, removing N 

