#!/bin/bash

### produce scripts to download and process public data with nfcore RNAseq pipeline: ###
### 1. FetchNGS  2. nfcore-rnaseq  3. differentialAbundance ###
### by Carlos Perez-Cervantes. Moskowitz lab. 2024-10-28 ###

# Parse command-line options
while getopts n:h arg; do
  case "$arg" in
    n) n="$OPTARG" ;;
    h)
      echo "Usage:"
      echo "   prep_nfcore_rnaseq_GEO2DApipeline.sh -h   Display this help message."
      echo "   -n: GEO ID"
      exit 0
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# Check if GEO ID is provided
if [ -z "$n" ]; then
  echo "Error: GEO ID not provided. Use -n to specify the GEO ID."
  exit 1
fi

mkdir -p /scratch/midway3/$USER/$n 

cd /scratch/midway3/$USER/$n 

echo -e "\n\nwork directory created at /scratch/midway3/$USER/$n\n saving outputs there now\n\n"


# Write prep script
cat << "EOM" > "fetchNGS_prep_$n.sh"
#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: prep_nfcore_rnaseq_GEO2DApipeline.sh [-r true|false]"
    exit 1
}

# Default value for resume parameter
resume=false

# Parse command-line options
while getopts "r:" opt; do
    case $opt in
        r )
            resume=$OPTARG
            ;;
        \? )
            usage
            ;;
    esac
done

module load python/anaconda-2022.05
module load singularity
source activate /project/imoskowitz/shared/software/nf-core
export TMPDIR="/scratch/midway3/$USER"
export NXF_TEMP="/scratch/midway3/$USER"
export NXF_SINGULARITY_CACHEDIR="/project/imoskowitz/shared/software/singularityImages"

echo  mystudy  > ids.csv
nf_command="nextflow run nf-core/fetchngs \
--input ids.csv \
--outdir . \
-profile singularity \
--nf_core_pipeline rnaseq \
--download_method aspera"

# Append -resume if the resume parameter is true
if [ "$resume" = "true" ]; then
    nf_command="$nf_command -resume"
fi

# Execute the nextflow command
eval $nf_command

echo -e "\n\nmodified samplesheet for DA saved as DA_samplesheet.csv\n\n"

csvcut -c sample,sample_alias,sample_title,sample_description samplesheet/samplesheet.csv > DA_samplesheet.csv

echo -e "\n\nPlease add condition and replicate columns to the DA sample sheet!\n\n"
echo -e "id,variable,reference,target,blocking\n$n,condition,control,treated," > contrast.csv
echo -e "\n\nDummy contrast file was created, please edit accordingly to match DA samplesheet\n\n"

EOM

cat << "EOM" > "nfcorernaseq_prep_$n.sh"
#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 [-r true|false]"
    exit 1
}

# Default value for resume parameter
resume=false

# Parse command-line options
while getopts "r:" opt; do
    case $opt in
        r )
            resume=$OPTARG
            ;;
        \? )
            usage
            ;;
    esac
done

module load python/anaconda-2022.05
module load singularity
source activate /project/imoskowitz/shared/software/nf-core
export TMPDIR="/scratch/midway3/$USER"
export NXF_TEMP="/scratch/midway3/$USER"
export NXF_SINGULARITY_CACHEDIR="/project/imoskowitz/shared/software/singularityImages"


nf_command="nextflow run \
/project/imoskowitz/shared/software/nf-core-rnaseq_3.13.2/3_13_2/ \
-work-dir temp \
-profile singularity \
--input samplesheet/samplesheet.csv \
--outdir nf_rnaseq_out  \
--genome mm10 \
--aligner star_salmon \
--gtf: /project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf \
--star_index: /project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/STARIndex \
--salmon_index: /project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/salmon_index \
--igenomes_base: /project/imoskowitz/shared/annotations/iGenomes/ \
--skip_rseqc true \
--skip_bbsplit true \
--skip_markduplicates true \
--skip_umi_extract true \
--skip_trimming true \
--skip_stringtie true \
--skip_dupradar true \
--skip_biotype_qc true \
--skip_deseq2_qc true \
--skip_bigwig true \
--skip_qualimap true \
--skip_fastqc true \
--skip_preseq true \
--skip_rseqc true \
-c /project/imoskowitz/shared/software/assets/test4.config"

# Append -resume if the resume parameter is true
if [ "$resume" = "true" ]; then
    nf_command="$nf_command -resume"
fi

# Execute the nextflow command
eval $nf_command

EOM

cat << "EOM" > "nfcoreDA_prep_$n.sh"
#!/bin/bash

module load python/anaconda-2022.05
module load singularity
source activate /project/imoskowitz/shared/software/nf-core
export TMPDIR="/scratch/midway3/$USER"
export NXF_TEMP="/scratch/midway3/$USER"
export NXF_SINGULARITY_CACHEDIR="/project/imoskowitz/shared/software/singularityImages"


nextflow run nf-core/differentialabundance \
--input DA_samplesheet.csv \
-work-dir /scratch/midway3/$USER/DAstudy \
--contrasts contrast.csv \
--study_name mystudy \
--matrix nf_rnaseq_out/star_salmon/salmon.merged.gene_counts_scaled.tsv \
--gtf /project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf \
--outdir DA \
-profile rnaseq,singularity \
--features_gtf_feature_type CDS \
--features_id_col gene_id  \
--deseq_vs_blind false \
--deseq2_alpha 0.05 \
--deseq2_vst_nsub 2000

EOM

sed -i "s/mystudy/$n/g" *.sh

echo -e "\n\nprep scripts created at /scratch/midway3/$USER/$n\n\n"



