#!/usr/bin/env python

##### helpful python function #########
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


def process_dataframe(dataframe):
    seqrunid = dataframe["Sequencing run ID"].drop(index=0).dropna()[2]
    expid = dataframe["Experiment ID"].drop(index=0).dropna()[2]
    dataframe = dataframe.drop(index=0)
    seqpath = "/project/imoskowitz/shared/sequencing.data"
    subdir=os.listdir(f"{seqpath}/{seqrunid}")[0]
    dataframe = dataframe[~dataframe.isin(['x']).any(axis=1)]
    if dataframe[["FASTQ.R2"]].isnull().values.all():
        dataframe["sample"] = dataframe["Experiment ID"] + "_" + dataframe["Sample ID"]
        dataframe["fastq_1"] = seqpath + "/" + dataframe['Sequencing run ID'] + "/" + subdir + "/" + dataframe['FASTQ.R1']
        dataframe = dataframe[['sample', 'fastq_1']]
        dataframe = dataframe[~dataframe.isin([np.nan]).any(axis=1)]
        dataframe = dataframe.assign(fastq_2="")
        dataframe = dataframe.assign(strandedness="auto")
    else:
        dataframe["sample"] = dataframe["Experiment ID"] + "_" + dataframe["Sample ID"]
        dataframe["fastq_1"] = seqpath + "/" + dataframe['Sequencing run ID'] + "/" + subdir + "/" + dataframe['FASTQ.R1']
        dataframe["fastq_2"] = seqpath + "/" + dataframe['Sequencing run ID'] + "/" + subdir + "/" + dataframe['FASTQ.R2']
        dataframe = dataframe[['sample', 'fastq_1', 'fastq_2']]
        dataframe = dataframe[~dataframe.isin([np.nan]).any(axis=1)]
        dataframe = dataframe.assign(strandedness="auto")
    return dataframe,seqrunid,expid


def update_yaml(yaml_file, input_path, output_dir, out_yaml_file):
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)
        data["input"] = input_path
        data["outdir"] = output_dir
        data["email"] = os.getlogin()+"@uchicago.edu"
    with open(out_yaml_file, "w") as f:
        yaml.dump(data, f, default_flow_style=False)


def slurm_writer(pbsDir, seed, processors, nodes,memory,partition, walltime, module, shell):
    pbsf=open(pbsDir+"/"+seed+'.sh','w')
    pbsf.write('#!/bin/bash' +'\n')
    pbsf.write('#SBATCH --job-name='+seed+ '\n')
    pbsf.write('#SBATCH -n '+processors+ '\n')
    pbsf.write('#SBATCH -N '+nodes+ '\n')
    pbsf.write('#SBATCH -t '+walltime+':00:00' + '\n')
    pbsf.write('#SBATCH --mem-per-cpu='+memory+'gb '+ '\n')
    pbsf.write('#SBATCH --mail-user='+os.getlogin()+'@uchicago.edu'+'\n')
    pbsf.write('#SBATCH --mail-type=end'+'\n')
    pbsf.write('#SBATCH --account=pi-imoskowitz' + '\n')
    pbsf.write('#SBATCH --partition='+partition+'\n')
    pbsf.write(module+' \n')
    pbsf.write(shell)
    pbsf.close()


def process_sample_sheet(sample_sheet):
    sample_sheet.columns = sample_sheet.columns.str.replace('Sample name.*', 'Sample name', regex=True)
    sample_sheet = sample_sheet[~sample_sheet.isin(['x']).any(axis=1)].iloc[2:]
    sample_sheet["sampleID"] = sample_sheet["Experiment ID"] + "_" + sample_sheet["Sample ID"]
    sample_sheet['replicate'] = sample_sheet['Sample name'].str.split('_').str[-1].str.lstrip('0')
    sample_sheet['sample'] = sample_sheet['Sample name'].str.rsplit('_', 1).str[0]
    sample_sheet = sample_sheet[['sample','sampleID','replicate']].dropna()
    sample_sheet['replicate'] = sample_sheet['replicate'].astype(int).astype(str)
    return sample_sheet


def process_atac_sheet(dataframe):
    seqrunid = dataframe["Sequencing run ID"].drop(index=0).dropna()[2]
    expid = dataframe["Experiment ID"].drop(index=0).dropna()[2]
    dataframe = dataframe.drop(index=0)
    seqpath = "/project/imoskowitz/shared/sequencing.data"
    subdir=os.listdir(f"{seqpath}/{seqrunid}")[0]
    dataframe = dataframe[~dataframe.isin(['x']).any(axis=1)]
    if dataframe[["FASTQ.R2"]].isnull().values.all():
        dataframe["sampleID"] = dataframe["Experiment ID"] + "_" + dataframe["Sample ID"]
        dataframe["fastq_1"] = seqpath + "/" + dataframe['Sequencing run ID'] + "/" + subdir + "/" + dataframe['FASTQ.R1']
        dataframe = dataframe[['sampleID', 'fastq_1']]
        dataframe = dataframe[~dataframe.isin([np.nan]).any(axis=1)]
        dataframe = dataframe.assign(fastq_2="")
    else:
        dataframe["sampleID"] = dataframe["Experiment ID"] + "_" + dataframe["Sample ID"]
        dataframe["fastq_1"] = seqpath + "/" + dataframe['Sequencing run ID'] + "/" + subdir + "/" + dataframe['FASTQ.R1']
        dataframe["fastq_2"] = seqpath + "/" + dataframe['Sequencing run ID'] + "/" + subdir + "/" + dataframe['FASTQ.R2']
        dataframe = dataframe[['sampleID', 'fastq_1', 'fastq_2']]
        dataframe = dataframe[~dataframe.isin([np.nan]).any(axis=1)]
    return dataframe,seqrunid,expid


def process_chip_sample_sheet_old(sample_sheet):
    sample_sheet.columns = sample_sheet.columns.str.replace('Sample name.*', 'Sample name', regex=True)
    sample_sheet = sample_sheet[~sample_sheet.isin(['x']).any(axis=1)].iloc[2:]
    sample_sheet["sampleID"] = sample_sheet["Experiment ID"] + "_" + sample_sheet["Sample ID"]
    sample_sheet['replicate'] = sample_sheet['Sample name'].str.split('_').str[-1]
    sample_sheet['sample'] = sample_sheet['Sample name'].str.rsplit('_', 1).str[0]
    sample_sheet['antibody'] = ''
    sample_sheet['control'] = ''
    sample_sheet['control_replicate'] = sample_sheet['replicate']
    sample_sheet = sample_sheet[['sample','sampleID','replicate','antibody','control','control_replicate']].dropna()
    mask_input = sample_sheet['sample'].str.contains('input', case=False, na=False)
    sample_sheet.loc[~mask_input, 'antibody'] = 'FLAG'
    sample_sheet.loc[~mask_input, 'control'] = sample_sheet.loc[mask_input, 'sample'].iloc[0]
    sample_sheet = sample_sheet[['sample','sampleID','replicate','antibody','control','control_replicate']].dropna()
    return sample_sheet

def process_chip_sample_sheet_220(sample_sheet):
    sample_sheet.columns = sample_sheet.columns.str.replace('Sample name.*', 'Sample name', regex=True)
    sample_sheet = sample_sheet[~sample_sheet.isin(['x']).any(axis=1)].iloc[2:]
    sample_sheet["sampleID"] = sample_sheet["Experiment ID"] + "_" + sample_sheet["Sample ID"]
    sample_sheet['replicate'] = sample_sheet['Sample name'].str.split('_').str[-1]
    sample_sheet['sample'] = sample_sheet['Sample name'].str.rsplit('_', 1).str[0]
    sample_sheet['antibody'] = ''
    sample_sheet['control'] = ''
    sample_sheet['control_replicate'] = sample_sheet['replicate']
    sample_sheet = sample_sheet[['sample','sampleID','replicate','antibody','control','control_replicate']].dropna()
    mask_input = sample_sheet['sample'].str.contains('input', case=False, na=False)
    sample_sheet.loc[~mask_input, 'antibody'] = 'FLAG'
    for idx in sample_sheet[~mask_input].index:
        sample_name = sample_sheet.loc[idx, 'sample']
        sample_rep = sample_sheet.loc[idx, 'replicate']
        expected_input_name = sample_name.replace('FLAG', 'input').replace('flag', 'input')
        matching_inputs = sample_sheet[
            mask_input & 
            sample_sheet['sample'].str.contains(expected_input_name.split('_')[0], case=False)
        ]
        if len(matching_inputs) > 1:
            exact_matches = matching_inputs[
                matching_inputs['sample'].str.lower() == expected_input_name.lower()
            ]
            if not exact_matches.empty:
                matching_inputs = exact_matches
            else:
                for part in expected_input_name.split('_')[1:]:
                    if part.lower() not in ['flag', 'input']:
                        better_matches = matching_inputs[
                            matching_inputs['sample'].str.contains(part, case=False)
                        ]
                        if not better_matches.empty:
                            matching_inputs = better_matches
        if not matching_inputs.empty:
            rep_match = matching_inputs[matching_inputs['replicate'] == sample_rep]
            if not rep_match.empty:
                sample_sheet.loc[idx, 'control'] = rep_match['sample'].iloc[0]
                sample_sheet.loc[idx, 'control_replicate'] = rep_match['replicate'].iloc[0]
            else:
                sample_sheet.loc[idx, 'control'] = matching_inputs['sample'].iloc[0]
                sample_sheet.loc[idx, 'control_replicate'] = matching_inputs['replicate'].iloc[0]
    sample_sheet = sample_sheet[['sample','sampleID','replicate','antibody','control','control_replicate']].dropna()
    return sample_sheet



def process_chip_sample_sheet_200(sample_sheet):
    sample_sheet.columns = sample_sheet.columns.str.replace('Sample name.*', 'Sample name', regex=True)
    sample_sheet = sample_sheet[~sample_sheet.isin(['x']).any(axis=1)].iloc[2:]
    sample_sheet["sampleID"] = sample_sheet["Experiment ID"] + "_" + sample_sheet["Sample ID"]
    sample_sheet['replicate'] = sample_sheet['Sample name'].str.split('_').str[-1]
    sample_sheet['sample'] = sample_sheet['Sample name'].str.rsplit('_', 1).str[0]
    sample_sheet['antibody'] = ''
    sample_sheet['control'] = ''
    sample_sheet['control_replicate'] = sample_sheet['replicate']
    sample_sheet = sample_sheet[['sample','sampleID','replicate','antibody','control','control_replicate']].dropna()
    mask_input = sample_sheet['sample'].str.contains('input', case=False, na=False)
    sample_sheet.loc[~mask_input, 'antibody'] = 'FLAG'
    for idx in sample_sheet[~mask_input].index:
        sample_name = sample_sheet.loc[idx, 'sample']
        sample_rep = sample_sheet.loc[idx, 'replicate']
        expected_input_name = sample_name.replace('FLAG', 'input').replace('flag', 'input')
        matching_inputs = sample_sheet[
            mask_input & 
            sample_sheet['sample'].str.contains(expected_input_name.split('_')[0], case=False)
        ]
        if len(matching_inputs) > 1:
            exact_matches = matching_inputs[
                matching_inputs['sample'].str.lower() == expected_input_name.lower()
            ]
            if not exact_matches.empty:
                matching_inputs = exact_matches
            else:
                for part in expected_input_name.split('_')[1:]:
                    if part.lower() not in ['flag', 'input']:
                        better_matches = matching_inputs[
                            matching_inputs['sample'].str.contains(part, case=False)
                        ]
                        if not better_matches.empty:
                            matching_inputs = better_matches
        if not matching_inputs.empty:
            rep_match = matching_inputs[matching_inputs['replicate'] == sample_rep]
            if not rep_match.empty:
                sample_sheet.loc[idx, 'control'] = rep_match['sample'].iloc[0]
                sample_sheet.loc[idx, 'control_replicate'] = rep_match['replicate'].iloc[0]
            else:
                sample_sheet.loc[idx, 'control'] = matching_inputs['sample'].iloc[0]
                sample_sheet.loc[idx, 'control_replicate'] = matching_inputs['replicate'].iloc[0]
    sample_sheet = sample_sheet[['sample','sampleID','replicate','antibody','control','control_replicate']].dropna()
    sample_sheet['sample'] = sample_sheet['sample']+'_REP'+sample_sheet['replicate']
    sample_sheet = sample_sheet[['sample','sampleID','antibody','control']].dropna()
    return sample_sheet