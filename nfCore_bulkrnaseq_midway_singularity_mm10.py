#!/usr/bin/env python

##### run nf core rnaseq#########
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

from carlos_nfcFunctions import process_dataframe, update_yaml, slurm_writer

###########################################################
##
##  define parser and parameters for user input
###########################################################

parser = argparse.ArgumentParser(add_help=True)

parser.add_argument("-s", "--msheet", type=str,nargs="+",
		help="metasheet input to read")
parser.add_argument("-o", "--out", type=str, default="/project/imoskowitz/shared/sequencing.processed",
    help="specify output destination, but you better not unless youre sure!")

args = parser.parse_args()
yaml_file = "/project/imoskowitz/shared/software/assets/nfcore_bulkRNAseq_mm10.yaml"


if len(args.msheet) < 2:
    dataframe = pd.read_excel(args.msheet[0])
    processed_dataframe,seqrunid,expid = process_dataframe(dataframe)
    if not os.path.exists(args.out+"/"+seqrunid):
        os.makedirs(args.out+"/"+seqrunid)
    processed_dataframe.to_csv(args.out+"/"+seqrunid+"/"+"sample.sheet.csv", sep=',', index=False)
    input_path = args.out+"/"+seqrunid+"/"+"sample.sheet.csv"
    output_dir = args.out+"/"+seqrunid
    out_yaml_file = args.out+"/"+seqrunid+'/rnaseq_params.yaml'
    update_yaml(yaml_file, input_path, output_dir, out_yaml_file)
    slurm_writer(pbsDir = args.out+"/"+seqrunid, 
    seed = "nfcore_bulk_rnaseq_run_"+date.today().strftime("%y%m%d"), 
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
    "/project/imoskowitz/shared/software/nf-core-rnaseq_3.13.2/3_13_2/ -work-dir "+'/scratch/midway3/'+os.getlogin()+"/working_"+seqrunid+ " \\\n"+
    "-profile singularity -params-file "+args.out+"/"+seqrunid+"/rnaseq_params.yaml -c /project/imoskowitz/shared/software/assets/test4.config \n")
    print("\n\nnfcore files have been sent to "+args.out+"/"+seqrunid+"\n\n")

else:
    _xlist = []
    for sheet in args.msheet:
        dataframe = pd.read_excel(sheet)
        dataframex,seqrunid,expid = process_dataframe(dataframe)
        _xlist.append(dataframex)
    processed_dataframe = pd.concat(_xlist)
    if not os.path.exists(args.out+"/"+expid):
        os.makedirs(args.out+"/"+expid)
    processed_dataframe.to_csv(args.out+"/"+expid+"/"+"sample.sheet.csv", sep=',', index=False)
    input_path = args.out+"/"+expid+"/"+"sample.sheet.csv"
    output_dir = args.out+"/"+expid
    out_yaml_file = args.out+"/"+expid+'/rnaseq_params.yaml'
    update_yaml(yaml_file, input_path, output_dir, out_yaml_file)
    slurm_writer(pbsDir = args.out+"/"+expid, 
    seed = "nfcore_bulk_rnaseq_run_"+date.today().strftime("%y%m%d"), 
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
    "/project/imoskowitz/shared/software/nf-core-rnaseq_3.13.2/3_13_2/ -work-dir "+'/scratch/midway3/'+os.getlogin()+"/working_"+expid+ " \\\n"+
    "-profile singularity -params-file "+args.out+"/"+expid+"/rnaseq_params.yaml -c /project/imoskowitz/shared/software/assets/test4.config \n")
    print("\n\nnfcore files have been sent to "+args.out+"/"+expid+"\n\n")
