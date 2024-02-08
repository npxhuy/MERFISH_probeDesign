# Name of file: loop_paintshop.sh

# Create by Hy on Nov 2023

# Description: This script is used to run PaintSHOP pipeline for multiple genomes in a folder

# Requirment: Python 3, Snakemake, Conda
# paintshop_snakemake env must be install in conda before running this script


source activate paintshop_snakemake
#To go to the folder containing the data, maybe can change to cd ../data/copy_extract_data but since it works with this, I will leave it like this.
# For the lunarc version, this does not require to go to the folder containing the data, cause we have flag to show the directory
cd ..
ls
cd data/extract_data
ls | while read folder; do
  
  cd $folder

  # Make config.yml file
  echo "assembly: ""'"$folder"'" > config.yml
  echo "genome_fasta: " "'"$folder".fna'" >> config.yml
  echo "annotation_file: " "'"$folder".gtf'" >> config.yml
  echo "blockparse_min_length: 30" >> config.yml
  echo "blockparse_max_length: 30" >> config.yml
  echo "blockparse_min_tm: 42" >> config.yml
  echo "blockparse_max_tm: 47" >> config.yml
  echo "model_temp: 42" >> config.yml


  # Run PaintSHOP pipeline
  # configure file paths
  CONFIG_FILE='config.yml'  
  SNAKE_FILE='../../../paintshop/PaintSHOP_pipeline/workflow/Snakefile'
  CONDA_ENVS='../../../paintshop/PaintSHOP_pipeline/shared_conda_envs'
  
  # run the pipeline


  cd ..
done


module load Anaconda3
source activate /home/npxhuy/.conda/envs/base1/envs/paintshop_snakemake