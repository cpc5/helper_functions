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

###########################################################
##
##  define parser and parameters for user input
###########################################################

parser = argparse.ArgumentParser(add_help=True)
	
def readable_dir(string):
	if not os.path.isdir(string):
		raise argparse.ArgumentTypeError("%s directory not found" % string) # allows a directory to be used as input
	return os.path.abspath(string)

parser.add_argument("-s", "--msheet", type=str,nargs="+",
		help="metasheet input to read")

parser.add_argument("-o", "--out", type=str, default="/project/imoskowitz/shared/sequencing.processed",
    help="specify output destination, but you better not unless youre sure!")

# parser.add_argument("-a", "--otherOps", type=str,default='--aligner "bowtie2" ',
#         help="other option to consider. ie --aligner 'bowite2' see manual at https://nf-co.re/chipseq/2.0.0/usage ")

args = parser.parse_args()

### to do change so that this intakes multiples sheet and combines them
if len(args.msheet) < 2:
    dataframe = pd.read_excel(args.msheet[0])

seqpath = "/project/imoskowitz/shared/sequencing.data"

seqrunid = dataframe["Sequencing run ID"].drop(index=0).dropna()[2]

if not os.path.exists(args.out+"/"+seqrunid):
    os.makedirs(args.out+"/"+seqrunid)

dataframe = dataframe.drop(index=0)
dataframe = dataframe[~dataframe.isin(['x']).any(axis=1)]

#dataframe = dataframe["Sample_name"] = dataframe["Experiment ID"] +"_"+dataframe["Sample ID"]

if dataframe[["FASTQ.R2"]].isnull().values.all():
    dataframe["sample"] = dataframe["Experiment ID"] +"_"+dataframe["Sample ID"]
    dataframe["fastq_1"] = seqpath+"/"+dataframe['Sequencing run ID']+"/FastX/"+dataframe['FASTQ.R1']
    dataframe = dataframe[['sample', 'fastq_1']]
    dataframe = dataframe[~dataframe.isin([np.nan]).any(axis=1)]
    dataframe = dataframe.assign(fastq_2 = "")
    dataframe = dataframe.assign(strandedness = "auto")
else:
    dataframe["sample"] = dataframe["Experiment ID"] +"_"+dataframe["Sample ID"]
    dataframe["fastq_1"] = seqpath+"/"+dataframe['Sequencing run ID']+"/FastX/"+dataframe['FASTQ.R1']
    dataframe["fastq_2"] = seqpath+"/"+dataframe['Sequencing run ID']+"/FastX/"+dataframe['FASTQ.R2']
    dataframe = dataframe[['sample', 'fastq_1', 'fastq_2']]
    dataframe = dataframe[~dataframe.isin([np.nan]).any(axis=1)]
    dataframe = dataframe.assign(strandedness = "auto")

dataframe.to_csv(args.out+"/"+seqrunid+"/"+"sample.sheet.csv", sep=',', index=False)


#####################
## edit YAML
#####################

def update_yaml(yaml_file, input_path, output_dir, out_yaml_file):
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)
        data["input"] = input_path
        data["outdir"] = output_dir
        data["email"] = os.getlogin()+"@uchicago.edu"
    with open(out_yaml_file, "w") as f:
        yaml.dump(data, f, default_flow_style=False)

yaml_file = "/project/imoskowitz/shared/software/assets/nfcore_bulkRNAseq_params_human.yaml"
input_path = args.out+"/"+seqrunid+"/"+"sample.sheet.csv"
output_dir = args.out+"/"+seqrunid
out_yaml_file = args.out+"/"+seqrunid+'/rnaseq_params.yaml'

update_yaml(yaml_file, input_path, output_dir, out_yaml_file)

#######################################################
##
##  create and define job dependencies
##
##########################################################

def slurm_writer(pbsDir, seed, processors, nodes,memory,partition, walltime, module, shell):
    pbsf=open(pbsDir+"/"+seed+'.sh','w')
    pbsf.write('#!/bin/bash' +'\n')
    pbsf.write('#SBATCH --job-name='+seed+ '\n')
    pbsf.write('#SBATCH -n '+nodes+ '\n')
    pbsf.write('#SBATCH -N '+processors+ '\n')
    pbsf.write('#SBATCH -t '+walltime+':00:00' + '\n')
    pbsf.write('#SBATCH --mem-per-cpu='+memory+'gb '+ '\n')
    pbsf.write('#SBATCH --mail-user='+os.getlogin()+'@uchicago.edu'+'\n')
    pbsf.write('#SBATCH --mail-type=end'+'\n')
    pbsf.write('#SBATCH --account=pi-imoskowitz' + '\n')
    pbsf.write('#SBATCH --partition='+partition+'\n')
    pbsf.write(module+' \n')
    pbsf.write(shell)
    pbsf.close()

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
"-profile singularity -params-file "+args.out+"/"+seqrunid+"/rnaseq_params.yaml -c /project/imoskowitz/shared/software/assets/test3.config \n")


print("\n\nnfcore files have been sent to "+args.out+"/"+seqrunid+"\n\n")

