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
  SNAKE_FILE='../../../04_tools/PaintSHOP_pipeline/workflow/Snakefile'
  CONDA_ENVS='../../../04_tools/PaintSHOP_pipeline/shared_conda_envs'
  
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