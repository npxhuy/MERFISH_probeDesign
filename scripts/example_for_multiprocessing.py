#!/bin/python
#SBATCH -A your_project
#SBATCH -p your_partition
#SBATCH -t your:time:limit
#SBATCH -N number_of_nodes
#SBATCH --ntasks-per-node = number_of_cores
#SBATCH --mail-user = name@yourmail.com
#SBATCH --mail-type = ALL
#SBATCH -J your_job_name
#SBATCH -o /directory/for/output/name.out
#SBATCH -e /directory/for/error/name.err
#SBATCH -D /directory/that/has/the/fasta/file/that/you/want/to/blast

# Sample code for multiprocessing on blastn, you have to modify a lot of things to make it work for your case like the directory paths, the number of files to split, the number of cores, the blastn path, the database path, etc.

import subprocess
import multiprocessing
import sys

# Define the path to the blastn and the human RNA database
blastn = "/path/to/blastn/ncbi-blast-2.15.0+/bin/blastn"
db = "/path/to/human/db/GRCh38_rna"

def run_blastn(file_name):
    command = f"{blastn} -db {db}  -query split_files/{file_name} -word_size 11 -ungapped > blast_result/{file_name}"
    subprocess.run(command, shell=True)

if __name__ == '__main__':
    # Split the file
    file_name = sys.argv[1]
    command =f"mkdir split_files; mkdir blast_result; split -d -l $(($(wc -l < {file_name}) / 39 )) {file_name} split_files/{file_name}"
    subprocess.run(command, shell=True)
    # List of file names
    files = [f"{file_name}{str(i).zfill(2)}" for i in range(40)]
    
    # Create a pool of worker processes
    with multiprocessing.Pool() as pool:
        pool.map(run_blastn, files)