#!/bin/bash

### RGT-hint footprinting analysis from ATAC-seq data ###
### Carlos Perez-Cervantes Moskowitz lab

# Parse command-line options
while getopts t:c:p:q:o:n:h arg; do
  case "$arg" in
    t) treat_bam="$OPTARG" ;;
    c) control_bam="$OPTARG" ;;
    p) treat_peaks="$OPTARG" ;;
    q) control_peaks="$OPTARG" ;;
    o) organism="$OPTARG" ;;
    n) ncpu="$OPTARG" ;;
    h)
      echo "Usage: rgt_hint_footprinting.sh -t <treatment_bam> -c <control_bam> -p <treatment_peaks> -q <control_peaks> -o <organism> -n <cpu_cores>"
      echo "   -t: Treatment BAM file from ATAC-seq alignment"
      echo "   -c: Control BAM file from ATAC-seq alignment (optional for differential analysis)"
      echo "   -p: Treatment peaks file in narrowPeak format"
      echo "   -q: Control peaks file in narrowPeak format (optional for differential analysis)"
      echo "   -o: Organism (mouse or human)"
      echo "   -n: Number of CPU cores for parallel processing (default: 1)"
      echo "   -h: Display this help message"
      exit 0
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# Set default value for ncpu if not provided
if [ -z "$ncpu" ]; then
  ncpu=1
fi

# Check if required parameters are provided
if [ -z "$treat_bam" ] || [ -z "$treat_peaks" ] || [ -z "$organism" ]; then
  echo "Error: missing required parameters"
  echo "Required: -t <treatment_bam> -p <treatment_peaks> -o <organism>"
  echo "Use -h for help"
  exit 1
fi

# Set organism-specific parameters
case "$organism" in
  mouse|Mouse|MOUSE)
    rgt_organism="mm10"
    genome_fasta="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
    genome_gtf="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
    genome_bed="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.bed"
    blacklist="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/mm10-blacklist.v2.bed"
    ;;
  human|Human|HUMAN)
    rgt_organism="hg38"
    genome_fasta="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
    genome_gtf="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
    genome_bed="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed"
    blacklist="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/hg38-blacklist.v2.bed"
    ;;
  *)
    echo "Error: organism must be 'mouse' or 'human'"
    exit 1
    ;;
esac

module load singularity

# Define the RGT container
RGT_CONTAINER="/project/imoskowitz/shared/software/singularityImages/singularity-images/depot.galaxyproject.org-singularity-rgt-1.0.2--py37he4a0461_0.img"

# Extract sample names from BAM files (remove path and extension)
treat_prefix=$(basename "$treat_bam" .bam)
if [ ! -z "$control_bam" ]; then
  control_prefix=$(basename "$control_bam" .bam)
fi

# Create output directories
mkdir -p ./match
mkdir -p ./tracks
mkdir -p ./footprints
mkdir -p ./differential

echo -e "\nRunning RGT-hint footprinting on treatment sample\n"

singularity exec --bind "$PWD" \
"$RGT_CONTAINER" \
rgt-hint footprinting --atac-seq --paired-end --organism="$rgt_organism" \
--output-location=./footprints --output-prefix="$treat_prefix" \
"$treat_bam" "$treat_peaks" || { echo "Error in RGT-hint footprinting step for treatment"; exit 1; }

# Run footprinting on control if provided
if [ ! -z "$control_bam" ] && [ ! -z "$control_peaks" ]; then
  echo -e "\nRunning RGT-hint footprinting on control sample\n"
  
  singularity exec --bind "$PWD" \
  "$RGT_CONTAINER" \
  rgt-hint footprinting --atac-seq --paired-end --organism="$rgt_organism" \
  --output-location=./footprints --output-prefix="$control_prefix" \
  "$control_bam" "$control_peaks" || { echo "Error in RGT-hint footprinting step for control"; exit 1; }
fi

echo -e "\nRunning RGT-hint tracks for treatment sample\n"

singularity exec --bind "$PWD" \
"$RGT_CONTAINER" \
rgt-hint tracks --bc --bigWig --organism="$rgt_organism" \
"$treat_bam" "$treat_peaks" --output-prefix="$treat_prefix"_BC \
--output-location=./tracks || { echo "Error in RGT-hint tracks step for treatment"; exit 1; }

# Run tracks on control if provided
if [ ! -z "$control_bam" ] && [ ! -z "$control_peaks" ]; then
  echo -e "\nRunning RGT-hint tracks for control sample\n"
  
  singularity exec --bind "$PWD" \
  "$RGT_CONTAINER" \
  rgt-hint tracks --bc --bigWig --organism="$rgt_organism" \
  "$control_bam" "$control_peaks" --output-prefix="$control_prefix"_BC \
  --output-location=./tracks || { echo "Error in RGT-hint tracks step for control"; exit 1; }
fi

# Run motif analysis if both treatment and control are provided
if [ ! -z "$control_bam" ] && [ ! -z "$control_peaks" ]; then
  echo -e "\nRunning RGT motif analysis\n"
  
  singularity exec --bind "$PWD" \
  "$RGT_CONTAINER" \
  rgt-motifanalysis matching --organism="$rgt_organism" \
  --input-files ./footprints/"$treat_prefix".bed ./footprints/"$control_prefix".bed \
  --output-location=./match || { echo "Error in RGT motif analysis step"; exit 1; }
  
  echo -e "\nRunning RGT-hint differential analysis\n"
  
  singularity exec --bind "$PWD" \
  "$RGT_CONTAINER" \
  rgt-hint differential --organism="$rgt_organism" --bc --nc "$ncpu" \
  --mpbs-files=./match/"$treat_prefix"_mpbs.bed,./match/"$control_prefix"_mpbs.bed \
  --reads-files="$treat_bam","$control_bam" \
  --conditions="$treat_prefix","$control_prefix" \
  --output-location=./differential/"$treat_prefix"_"$control_prefix" || { echo "Error in RGT-hint differential analysis step"; exit 1; }
fi

echo -e "\nRGT-hint analysis completed successfully\n"
echo -e "Results are available in the following directories:\n"
echo -e "  - Footprints: ./footprints/"
echo -e "  - Tracks: ./tracks/"
if [ ! -z "$control_bam" ] && [ ! -z "$control_peaks" ]; then
  echo -e "  - Motif matching: ./match/"
  echo -e "  - Differential analysis: ./differential/"
fi 