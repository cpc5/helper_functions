#!/usr/bin/env python

##### WRITE PBS SCRIPT AND CALL STAR #########
#	This python script will write a .pbs for every fastq file      #
#	listed in a metasheet. 
#	By:Carlos PerezCervantes. Moskowitz lab # 
#####################################################

import os, datetime, time, argparse, importlib
import pandas as pd
from glob import glob
from  builtins import any as b_any


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

parser.add_argument("-d", "--fastqDir", type=str,nargs="+",
        help="directory of Fastq files")

parser.add_argument("-o", "--out", type=str, default="/scratch/midway3/"+os.getlogin()+"/temp",
		help="specify output destination")
				  
parser.add_argument("-p", "--processors", type=str, default='04',
    help="specify number of processors for computation node, ..02, 08")

parser.add_argument("-n", "--nodes", type=str, default='1',
    help="specify number of node, ..02, 08")

parser.add_argument("-w", "--walltime", type=str, default='08', 
    help="specify running time in hours, i.e 02, 4, 16")

parser.add_argument("-m", "--memory", type=str, default='64', 
    help="specify memory use in GB, i.e 02, 4, 16 ")

parser.add_argument("-l", "--libraryStrat", type=str, default='single', 
    help="specify 'single' or 'paired' end reads")

parser.add_argument("-f", "--partition", type=str, default='caslake', 
    help="specify partition to use; see [rcchelp]")

parser.add_argument("-g", "--genome", type=str,default='mm10',
        help="specify genome: mm10 or hg38")


args = parser.parse_args()
args.out = os.path.abspath(args.out)
#data = pd.read_excel(args.msheet)



##############################################################################


def pbs_writer(pbsDir, seed, processors, nodes,memory,partition, walltime, module, shell):
    pbsf=open(pbsDir+seed+'.sh','w')
    pbsf.write('#!/bin/bash' +'\n')
    pbsf.write('#SBATCH --job-name=' +seed+ '\n')
    pbsf.write('#SBATCH -N '+nodes+ '\n')
    pbsf.write('#SBATCH -n '+processors+ '\n')
    pbsf.write('#SBATCH -t '+walltime+':00:00' + '\n')
    pbsf.write('#SBATCH --mem-per-cpu='+args.memory+'gb '+ '\n')
    pbsf.write('#SBATCH --mail-user='+os.getlogin()+'@uchicago.edu'+'\n')
    pbsf.write('#SBATCH --mail-type=end'+'\n')
    pbsf.write('#SBATCH --account=pi-imoskowitz' + '\n')
    pbsf.write('#SBATCH --partition='+partition+'\n')
    pbsf.write(module+' \n')
    pbsf.write('cd $PBS_O_WORKDIR'+'\n\n')
    pbsf.write(shell)
    pbsf.close()

########################################################################

##########################################3
##
##   create log and pbs directory
##
###########################################

if not os.path.exists(args.out):
    os.makedirs(args.out)

jobsDir = args.out+'/log/'

if not os.path.exists(jobsDir):
    os.makedirs(jobsDir)

###########################################
##
##  Define genome locations from input
###########################################

# species = data["Species"][0]

# if species == 'mouse':

#     genome_index =' /gpfs/data/referenceFiles/Mus_musculus/STARgenome/UCSC-mm10-50bp '
#     genome =' /gpfs/data/referenceFiles/Mus_musculus/STARgenome/UCSC-mm10-50bp/Genome '
#     gtf = '/gpfs/data/referenceFiles/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf '
#     genomeSize = '2652783500'

# elif args.genome == 'human':

#     genome = '/gpfs/data/referenceFiles/Homo_sapiens/STARgenome/GRCh38.primary_Gencode27_50bp/Genome '
#     genome_index =' /gpfs/data/referenceFiles/Homo_sapiens/STARgenome/GRCh38.primary_Gencode27_50bp/ '
#     gtf = ' /gpfs/data/referenceFiles/Homo_sapiens/UCSC/hg38/Annotation/knownGene_hg38.gtf '
#     genomeSize = '2913022398'
# else:
#     print("no species column found")


if args.genome=='mm10':

    genome_index =' /project/imoskowitz/shared/annotations/STARgenome/GRCm38.p6_vM23_100bp/ '
    genome =' /project/imoskowitz/shared/annotations/STARgenome/GRCm38.p6_vM23_100bp/Genome '
    gtf= '/project/imoskowitz/shared/annotations/STARgenome/GRCm38.p6_vM23_100bp/gencode.vM23.annotation.gtf '
    genomeSize = '2652783500'
elif args.genome == 'hg38':

    genome = '/project/imoskowitz/shared/annotations/STARgenome/GRCh38_v24_100bp/Genome'
    genome_index ='/project/imoskowitz/shared/annotations/STARgenome/GRCh38_v24_100bp/STARIndex' 
    gtf='/project/imoskowitz/shared/annotations/STARgenome/GRCh38_v24_100bp/gencode.v24.primary_protein_coding.gtf'
    genomeSize = '2913022398'
else:
    print("specify genome")


##########################################################
##
##  create PBS and define job dependencies
##
##########################################################
# muhFiles = data["Directory path"]+"/Fastq/"+data["FASTQ.R1"]
# outputNames = data["Sample name"]
# outputNameList = outputNames.tolist()
# FilesList = muhFiles.tolist()

FilesList = []
for dir in args.fastqDir:
    for fastq in glob(os.path.join(os.path.abspath(dir), '*[!Undetermined].fastq.gz*')):
        FilesList.append(fastq)

#FilesList= glob(os.path.join(args.fastqDir, '*.fastq.gz'))
thingtemp=[os.path.basename(i) for i in FilesList]
outputNameList=list(set(x.replace(".fastq.gz", '').replace(".fastq.gz", '') for x in thingtemp))

myDir1=glob(os.path.join(os.path.abspath(args.fastqDir[0]), '*[!Undetermined].fastq.gz*'))
sampleName = []
for fastq in myDir1:
    name = os.path.basename(fastq).replace(".fastq.gz", '').strip().split("_")[0]
    sampleName.append(name)

for fs in FilesList:
    for xm in sampleName:
        #f=fs.strip().split('/')[-1].strip().split('.')[0]
        if not os.path.exists(args.out+'/'+xm):
            os.makedirs(args.out+'/'+xm)
        if args.libraryStrat == "paired":
            sets = [x for x in FilesList if xm in x]
            mate1reads = [h for h in sets if any(mate in h for mate in mater1)]
            mate2reads = [h for h in sets if any(mate in h for mate in mater2)]
            readFilesIn=(",".join(mate1reads)+ " \\\n"+",".join(mate2reads))   # add mate
        else:
            readFilesIn=(",".join([x for x in FilesList if xm in x]))
        pbs_writer(pbsDir = jobsDir, 
        seed = xm,
        nodes = args.nodes,
        walltime = str(args.walltime),
        partition = args.partition, 
        processors = str(args.processors),
        memory=str(args.memory),
        module='module load python/anaconda-2022.05\n'+
        'source activate /project/imoskowitz/shared/software/envs/legacyRNApipe\n\n\n'+
        'TMPDIR=/scratch/midway3/'+os.getlogin()+'\nexport TMPDIR\n', 
        shell = "STAR --runThreadN "+args.processors+" --genomeDir "+genome_index+" "+ 
            "\\\n --readFilesIn "+readFilesIn+" \\\n --readFilesCommand zcat "+
            "\\\n --sjdbGTFfile  "+gtf+" "+
            "\\\n --quantMode GeneCounts "+  
            "\\\n --outFileNamePrefix "+args.out+"/"+xm+"/  --outSAMtype BAM SortedByCoordinate\n\n"
            "samtools index -@ "+args.processors+" "+args.out+"/"+xm+"/Aligned.sortedByCoord.out.bam\n\n"
            )




