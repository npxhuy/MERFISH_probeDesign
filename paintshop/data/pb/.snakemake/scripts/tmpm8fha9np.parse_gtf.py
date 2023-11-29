
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/Users/hy/miniconda3/envs/paintshop_snakemake/lib/python3.8/site-packages', '/Users/hy/Documents/GitHub/thesis/paintshop/PaintSHOP_pipeline/workflow/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x959\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c\x06pb.gtf\x94\x8cApipeline_output/01_reference_files/01_chrom_names/chrom_names.txt\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x03gtf\x94K\x00N\x86\x94\x8c\x0bchrom_names\x94K\x01N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x15\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x1b)}\x94\x8c\x05_name\x94h\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh\x0fh\nh\x11h\x0bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c=pipeline_output/01_reference_files/05_annotation_files/pb.bed\x94\x8c=pipeline_output/02_intermediate_files/01_chrom_gtf_dataframes\x94e}\x94(h\r}\x94(\x8c\x03bed\x94K\x00N\x86\x94\x8c\x0cchrom_df_dir\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh-h)h/h*ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8c\x0330G\x94\x8c\x053:0:0\x94e}\x94(h\r}\x94(\x8c\x05mfree\x94K\x00N\x86\x94\x8c\x04h_rt\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bhBh>hDh?ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bheK\x01hgK\x01ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x08assembly\x94\x8c\x02pb\x94\x8c\x0cgenome_fasta\x94\x8c\x06pb.fna\x94\x8c\x0fannotation_file\x94\x8c\x06pb.gtf\x94u\x8c\x04rule\x94\x8c\tparse_gtf\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cO/Users/hy/Documents/GitHub/thesis/paintshop/PaintSHOP_pipeline/workflow/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/hy/Documents/GitHub/thesis/paintshop/PaintSHOP_pipeline/workflow/scripts/parse_gtf.py';
######## snakemake preamble end #########

import os
import pickle

import numpy as np
import pandas as pd

# configure file paths
CHROM_NAMES = snakemake.input.chrom_names
GTF_PATH = snakemake.input.gtf
BED_PATH = snakemake.output.bed
CHROM_DF_DIR = snakemake.output.chrom_df_dir

# ensure target directory exists
os.makedirs(CHROM_DF_DIR, exist_ok=True)

def main():

    # load the input annotation file
    df = load_gtf(GTF_PATH)

    # filter annotations to only exons on canonical chromosomes
    df = parse_and_filter(df)

    # split gtf dataframe into individual chromosome files
    chroms = df['seqid'].unique()
    for chrom in chroms:

        # filter gtf to only the current chromosome
        chrom_df = df[df['seqid'] == chrom]

        # save filtered gtf data to disk for this chromosome
        df_path = os.path.join(CHROM_DF_DIR, f'{chrom}_filtered_gtf.dat')
        pickle.dump(chrom_df, open(df_path, 'wb'))

    # save exon-resolved bed file
    print('Saving isoform-resolved BED annotation file...')
    bed_df = df[[
        'seqid', 
        'start', 
        'end', 
        'transcript_id', 
        'score', 
        'strand',
        'transcript_version',
        'gene_id',
    ]]
    bed_df.to_csv(BED_PATH, sep='\t', index=False, header=None)


def load_gtf(file_path):

    # read input gtf file and rename columns
    df = pd.read_csv(GTF_PATH, sep='\t', header=None)
    df.rename(columns={
        0: 'seqid',
        1: 'source',
        2: 'type',
        3: 'start',
        4: 'end',
        5: 'score',
        6: 'strand',
        7: 'phase',
        8: 'attributes',
    }, inplace=True)

    # success
    return(df)


def parse_and_filter(df):

    # load canonical chromosome names
    chrom_names = get_chrom_names()

    # filter annotations to canonical chromosomes only
    df = df[df['seqid'].isin(chrom_names)]

    # filter annotations to exon records only
    df = df[df['type'] == 'exon']

    # parse attributes string into a dict
    attr_data = df['attributes'].apply(parse_attributes)

    # extract attributes as new columns
    df['gene_id'] = attr_data.apply(lambda x: x['gene_id'])
    df['transcript_version'] = attr_data.apply(lambda x: x['transcript_id'])

    df['transcript_id'] = df['transcript_version'].apply(lambda x: x.split('.')[0])

    # df['exon_id'] = attr_data.apply(lambda x: x['exon_id'])

    # remove unnecessary column
    df.drop(['attributes'], axis=1, inplace=True)

    # success
    return(df)


def get_chrom_names():

    # load chromosome names from filtered fasta
    with open(CHROM_NAMES, 'r') as infile:
        text = infile.read().strip(' \n')
    chrom_names = text.split('\n')

    # success
    return(chrom_names)


def parse_attributes(data):

    # parse attributes
    split_data = data.replace('"', '').strip(' ;').split(';')
    attr_data = {}
    for field in split_data:
        split_field = field.strip(' ').split(' ')
        key = split_field[0]
        val = ' '.join(split_field[1:])
        attr_data[key] = val

    # success        
    return(attr_data)


if __name__ == '__main__':
    main()
