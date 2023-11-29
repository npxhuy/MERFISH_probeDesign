# Name of file: loop_paintshop.py

# Create by Hy on Nov 2023

# Description: This script is used to run PaintSHOP pipeline for multiple genomes in a folder

# Requirment: Python 3, Snakemake, Conda
# paintshop_snakemake env must be install in conda before running this script


source activate paintshop_snakemake
ls
cd data
ls | while read folder; do
  
  cd $folder

  # Make config.yml file
  echo "assembly: ""'"$folder"'" > config.yml
  echo "genome_fasta: " "'"$folder".fna'" >> config.yml
  echo "annotation_file: " "'"$folder".gtf'" >> config.yml
  
  # Run PaintSHOP pipeline
  # configure file paths
  CONFIG_FILE='config.yml'  
  SNAKE_FILE='../../PaintSHOP_pipeline/workflow/Snakefile'
  CONDA_ENVS='../../PaintSHOP_pipeline/shared_conda_envs'
  
  # run the pipeline
  snakemake --configfile config.yml --snakefile $SNAKE_FILE --use-conda --conda-prefix $CONDA_ENVS --cores --restart-times 3
  
  # export PDF and svg visualizations of the DAG structure of pipeline steps
  echo -e "Exporting pipeline DAG to svg and pdf..."
  snakemake --configfile config.yml --snakefile $SNAKE_FILE --dag > dag.dot
  dot -Tpdf dag.dot > pipeline_output/pipeline.pdf
  dot -Tsvg dag.dot > pipeline_output/pipeline.svg
  rm dag.dot
  
  echo -e "Generating pipeline HTML report..."
  snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --report pipeline_output/report.html --conda-frontend mamba

  cd ..
done
