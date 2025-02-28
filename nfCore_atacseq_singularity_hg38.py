#!/usr/bin/env python

##### run nf core atacseq#########
#	By:Carlos PerezCervantes. Moskowitz lab #
#####################################################
## module load python/anaconda-2023.09

import os, time, argparse, importlib
from datetime import date
import pandas as pd
from glob import glob
from  builtins import any as b_any
import numpy as np
import yaml

from carlos_nfcFunctions import process_sample_sheet, update_yaml, slurm_writer, process_atac_sheet

parser = argparse.ArgumentParser(add_help=True)

parser.add_argument("-s", "--msheet", type=str,nargs="+",
		help="metasheet input to read")
parser.add_argument("-a", "--ssheet", type=str,
		help="sample sheet input to read")
parser.add_argument("-o", "--out", type=str, default="/project/imoskowitz/shared/sequencing.processed",
    help="specify output destination, but you better not unless youre sure!")

args = parser.parse_args()
yaml_file = "/project/imoskowitz/shared/software/assets/nfcore_atacseq_params_hg38.yaml"

sample_sheet = pd.read_excel(args.ssheet)
sample_sheet = process_sample_sheet(sample_sheet)


if len(args.msheet) < 2:
    dataframe = pd.read_excel(args.msheet[0])
    processed_dataframe,seqrunid,expid = process_atac_sheet(dataframe)
    if not os.path.exists(args.out+"/"+seqrunid):
        os.makedirs(args.out+"/"+seqrunid)
    merged_df = pd.merge(processed_dataframe, sample_sheet, on='sampleID', how='left')
    if 'fastq_2' in merged_df.columns:
        merged_df = merged_df[['sample', 'fastq_1', 'fastq_2', 'replicate','sampleID']]
    else:
        merged_df = merged_df[['sample', 'fastq_1', 'replicate','sampleID']]
    merged_df.to_csv(args.out + "/" + seqrunid + "/" + "sample.sheet.csv", sep=',', index=False)
    input_path = args.out+"/"+seqrunid+"/"+"sample.sheet.csv"
    output_dir = args.out+"/"+seqrunid
    out_yaml_file = args.out+"/"+seqrunid+'/atacseq_params.yaml'
    update_yaml(yaml_file, input_path, output_dir, out_yaml_file)
    slurm_writer(pbsDir = args.out+"/"+seqrunid,
    seed = "nfcore_bulk_atacseq_run_"+date.today().strftime("%y%m%d"),
    nodes = str(1),
    walltime = str(24),
    processors = str(4),
    memory=str(64),
    partition="caslake",
    module='module load python/anaconda-2022.05\n'+
    "source activate /project/imoskowitz/shared/software/nf-core\n"+
    'module load singularity\n\n'+
    'export TMPDIR="/scratch/midway3/'+os.getlogin()+'"\n\n'+
    'export NXF_TEMP="/scratch/midway3/'+os.getlogin()+'"\n\n'+
    'export NXF_SINGULARITY_CACHEDIR="/project/imoskowitz/shared/software/singularityImages"',
    shell = "nextflow run \\\n"+
    "/project/imoskowitz/shared/software/nf-core-atacseq_2.1.2/2_1_2 -work-dir "+'/scratch/midway3/'+os.getlogin()+"/working_"+seqrunid+ " \\\n"+
    "-profile singularity -params-file "+args.out+"/"+seqrunid+"/atacseq_params.yaml -c /project/imoskowitz/shared/software/assets/test4.config \n")
    print("\n\nnfcore files have been sent to "+args.out+"/"+seqrunid+"\n\n")

else:
    _xlist = []
    for sheet in args.msheet:
        dataframe = pd.read_excel(sheet)
        dataframex,seqrunid,expid = process_atac_sheet(dataframe)
        _xlist.append(dataframex)
    processed_dataframe = pd.concat(_xlist)
    if not os.path.exists(args.out+"/"+expid):
        os.makedirs(args.out+"/"+expid)
    merged_df = pd.merge(processed_dataframe, sample_sheet, on='sampleID', how='left')
    if 'fastq_2' in merged_df.columns:
        merged_df = merged_df[['sample', 'fastq_1', 'fastq_2', 'replicate','sampleID']]
    else:
        merged_df = merged_df[['sample', 'fastq_1', 'replicate','sampleID']]
    merged_df.to_csv(args.out + "/" + expid + "/" + "sample.sheet.csv", sep=',', index=False)
    input_path = args.out+"/"+expid+"/"+"sample.sheet.csv"
    output_dir = args.out+"/"+expid
    out_yaml_file = args.out+"/"+expid+'/atacseq_params.yaml'
    update_yaml(yaml_file, input_path, output_dir, out_yaml_file)
    slurm_writer(pbsDir = args.out+"/"+expid,
    seed = "nfcore_bulk_atacseq_run_"+date.today().strftime("%y%m%d"),
    nodes = str(1),
    walltime = str(24),
    processors = str(4),
    memory=str(64),
    partition="caslake",
    module='module load python/anaconda-2022.05\n'+
    "source activate /project/imoskowitz/shared/software/nf-core\n"+
    'module load singularity\n\n'+
    'export TMPDIR="/scratch/midway3/'+os.getlogin()+'"\n\n'+
    'export NXF_TEMP="/scratch/midway3/'+os.getlogin()+'"\n\n'+
    'export NXF_SINGULARITY_CACHEDIR="/project/imoskowitz/shared/software/singularityImages"',
    shell = "nextflow run \\\n"+
    "/project/imoskowitz/shared/software/nf-core-atacseq_2.1.2/2_1_2 -work-dir "+'/scratch/midway3/'+os.getlogin()+"/working_"+expid+ " \\\n"+
    "-profile singularity -params-file "+args.out+"/"+expid+"/atacseq_params.yaml -c /project/imoskowitz/shared/software/assets/test4.config \n")
    print("\n\nnfcore files have been sent to "+args.out+"/"+expid+"\n\n")
