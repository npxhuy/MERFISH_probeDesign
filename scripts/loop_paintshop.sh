ls | while read folder; do
    cd $folder
    
    echo "assembly: ""'"$folder"'" > config.yml
    echo "genome_fasta: " "'"$folder".fna'" >> config.yml
    echo "annotation_file: " "'"$folder".gtf'" >> config.yml
    

    # configure file paths
    CONFIG_FILE='config.yml'
    SNAKE_FILE='../PaintSHOP_pipeline/workflow/Snakefile'
    CONDA_ENVS='../PaintSHOP_pipeline/shared_conda_envs'

    # activate conda environment
    source activate paintshop_snakemake

    # run the pipeline
    snakemake --configfile config.yml --snakefile $SNAKE_FILE \
        --use-conda --conda-prefix $CONDA_ENVS --cores \
        --restart-times 3

    # export PDF and svg visualizations of the DAG structure of pipeline steps
    echo -e "Exporting pipeline DAG to svg and pdf..."
    snakemake --configfile config.yml --snakefile $SNAKE_FILE --dag > dag.dot
    dot -Tpdf dag.dot > pipeline_output/pipeline.pdf
    dot -Tsvg dag.dot > pipeline_output/pipeline.svg
    rm dag.dot

    echo -e "Generating pipeline HTML report..."
    snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --report pipeline_output/report.html --conda-frontend mamba

    # success
    echo -e "\nDONE!\n"

    cd ..
    done

