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
shutil.copy('../tools/PaintSHOP_pipeline/example_run/config.yml', new_directory)
shutil.copy('../tools/PaintSHOP_pipeline/example_run/run_pipeline.sh', new_directory)

print(os.listdir(new_directory))
print(os.listdir(new_directory+'/data/'))