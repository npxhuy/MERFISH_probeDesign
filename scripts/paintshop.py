'''
paintshop.py

Description: 

User-defined function: 
Non-standard modules: os, errno, sys, Path, shutil

Procedure:
    
Input: 
Output: 

Usage: 

Version: 1.00
Date 2022-11-01
Name: Pham Xuan Huy Nguyen

'''





import os
import errno
#import sys
#from pathlib import Path
import shutil


assembly = 'SA' #nickname for this assembly
genome = '../data/short_testdata/sa.fna' # Path to genome fasta
annotation = '../data/short_testdata/sa.gtf' # Path to annotation gtf file
temperature = "37" # PaintSHOP supports model temperatures of 37, 42, 47, 52, and 60 C
new_directory = assembly + "_" + temperature # The name of the directory where this run will go

# Create a PaintSHOP directory for your run, then copy genome data
try:
  os.mkdir(new_directory)
  os.mkdir(new_directory+'/data')
  os.mkdir(new_directory+'/data/pipeline_output')
except OSError as e:
  if e.errno == errno.EEXIST:
    print('file exists')

shutil.copy(genome, new_directory+'/data/') 
shutil.copy(annotation, new_directory+'/data/')
shutil.copy('../tools/PaintSHOP_pipeline/example_run/config.yml', new_directory) # Unnecessary
shutil.copy('../tools/PaintSHOP_pipeline/example_run/run_pipeline.sh', new_directory) #Unneccessary

#print(os.listdir(new_directory))
#print(os.listdir(new_directory+'/data/'))

# This function creats config.yml and run_pipeline.sh files for a PaintSHOP run
# These simple files specify where to find files and how to run the pipeline

def GenerateConfigFiles(directory_path, assembly, genome_fasta, annotation_file, temp):
  # Write the config.yml for this temperature
  try: 
    os.remove(directory_path+'/config.yml')
    print("previous config.yml file deleted! ... New config.yml file created")
  except OSError as e:
    if e.errno == errno.EEXIST:
        print('New config.yml file created!')
  f = open(directory_path+'/config.yml', "w")
  f.write("assembly: " + "'" + assembly + "'")
  f.write('\n')
  f.write("genome_fasta: " + "'" + genome_fasta + "'")
  f.write('\n')
  f.write("annotation_file: " + "'" + annotation_file + "'")
  f.write('\n')
  f.write("model_temp: " +  "'" + temp + "'")
  f.close()

  # Write the run_pipeline.sh for this example
  try: 
    os.remove(directory_path+'/run_pipeline.sh')
    print("previous run_pipeline.sh file deleted! ... New run_pipeline.sh file created")
  except OSError as e:
    if e.errno == errno.EEXIST:
        print('New run_pipeline.sh file created!')
  f = open(directory_path+'/run_pipeline.sh', "w")


  f.write("""# configure file paths
  CONFIG_FILE='config.yml'
  SNAKE_FILE='../../tools/PaintSHOP_pipeline/workflow/Snakefile'
  CONDA_ENVS='../../tools/PaintSHOP_pipeline/shared_conda_envs'

  # activate conda environment
  source activate paintshop_snakemake

  # run the pipeline
  snakemake --configfile config.yml --snakefile $SNAKE_FILE \
    --use-conda --conda-prefix $CONDA_ENVS --cores 2\
    --restart-times 3

  # export PDF and svg visualizations of the DAG structure of pipeline steps
  echo -e "Exporting pipeline DAG to svg and pdf..."
  snakemake --configfile config.yml --snakefile $SNAKE_FILE --dag > dag.dot

  dot -Tpdf dag.dot > pipeline_output/pipeline.pdf
  dot -Tsvg dag.dot > pipeline_output/pipeline.svg
  rm dag.dot

  echo -e "Generating pipeline HTML report..."
  snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --report pipeline_output/report.html
  
  # success
  echo -e "\nDONE!\n"
  """)

  f.close()
  os.chmod(directory_path+'/run_pipeline.sh', 0o755)

GenerateConfigFiles(new_directory, assembly, genome, annotation, temperature)