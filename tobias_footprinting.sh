#!/bin/bash

### TOBIAS footprinting analysis from ATAC-seq data ###
### Carlos Perez-Cervantes Moskowitz lab

# Parse command-line options
while getopts t:c:p:q:o:m:n:h arg; do
  case "$arg" in
    t) treat_bam="$OPTARG" ;;
    c) control_bam="$OPTARG" ;;
    p) treat_peaks="$OPTARG" ;;
    q) control_peaks="$OPTARG" ;;
    o) organism="$OPTARG" ;;
    m) motifs="$OPTARG" ;;
    n) ncpu="$OPTARG" ;;
    h)
      echo "Usage: tobias_footprinting.sh -t <treatment_bam> -c <control_bam> -p <treatment_peaks> -q <control_peaks> -o <organism> -m <motif_file> -n <cpu_cores>"
      echo "   -t: Treatment BAM file from ATAC-seq alignment"
      echo "   -c: Control BAM file from ATAC-seq alignment (optional)"
      echo "   -p: Treatment peaks file in BED format"
      echo "   -q: Control peaks file in BED format (optional)"
      echo "   -o: Organism (mouse or human)"
      echo "   -m: Motif file (optional, default: JASPAR2020 vertebrates)"
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
    genome_fasta="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
    genome_gtf="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
    genome_bed="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.bed"
    blacklist="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/mm10-blacklist.v2.bed"
    chrom_sizes="/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Sequence/mm10.chrom.sizes"
    ;;
  human|Human|HUMAN)
    genome_fasta="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
    genome_gtf="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
    genome_bed="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed"
    blacklist="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/hg38-blacklist.v2.bed"
    chrom_sizes="/project/imoskowitz/shared/annotations/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/hg38.chrom.sizes"
    ;;
  *)
    echo "Error: organism must be 'mouse' or 'human'"
    exit 1
    ;;
esac

# Set default motif file if not provided
if [ -z "$motifs" ]; then
  motifs="/project/imoskowitz/shared/annotations/motifs/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
fi

module load singularity

# Define the TOBIAS container
TOBIAS_CONTAINER="/project/imoskowitz/shared/software/singularityImages/singularity-images/depot.galaxyproject.org-singularity-tobias-0.17.1--py39hff726c5_0.img"

# Extract sample names from BAM files (remove path and extension)
treat_prefix=$(basename "$treat_bam" .bam)
if [ ! -z "$control_bam" ]; then
  control_prefix=$(basename "$control_bam" .bam)
fi

# Create output directory if it doesn't exist
mkdir -p "${treat_prefix}_results"

echo -e "\nRunning TOBIAS ATACorrect on treatment sample to correct Tn5 bias\n"

singularity exec --bind "$PWD" \
"$TOBIAS_CONTAINER" \
TOBIAS ATACorrect \
--bam "$treat_bam" \
--genome "$genome_fasta" \
--peaks "$treat_peaks" \
--blacklist "$blacklist" \
--outdir "${treat_prefix}_results" \
--prefix "${treat_prefix}" || { echo "Error in TOBIAS ATACorrect step for treatment"; exit 1; }

# Run ATACorrect on control if provided
if [ ! -z "$control_bam" ] && [ ! -z "$control_peaks" ]; then
  echo -e "\nRunning TOBIAS ATACorrect on control sample to correct Tn5 bias\n"
  
  singularity exec --bind "$PWD" \
  "$TOBIAS_CONTAINER" \
  TOBIAS ATACorrect \
  --bam "$control_bam" \
  --genome "$genome_fasta" \
  --peaks "$control_peaks" \
  --blacklist "$blacklist" \
  --outdir "${treat_prefix}_results" \
  --prefix "${control_prefix}" || { echo "Error in TOBIAS ATACorrect step for control"; exit 1; }
fi

echo -e "\nRunning TOBIAS FootprintScores on treatment sample to calculate footprint scores\n"

singularity exec --bind "$PWD" \
"$TOBIAS_CONTAINER" \
TOBIAS FootprintScores \
--signal "${treat_prefix}_results/${treat_prefix}_corrected.bw" \
--peaks "$treat_peaks" \
--genome "$genome_fasta" \
--output "${treat_prefix}_results/${treat_prefix}_footprints.bw" || { echo "Error in TOBIAS FootprintScores step for treatment"; exit 1; }

# Run FootprintScores on control if provided
if [ ! -z "$control_bam" ] && [ ! -z "$control_peaks" ]; then
  echo -e "\nRunning TOBIAS FootprintScores on control sample to calculate footprint scores\n"
  
  singularity exec --bind "$PWD" \
  "$TOBIAS_CONTAINER" \
  TOBIAS FootprintScores \
  --signal "${treat_prefix}_results/${control_prefix}_corrected.bw" \
  --peaks "$control_peaks" \
  --genome "$genome_fasta" \
  --output "${treat_prefix}_results/${control_prefix}_footprints.bw" || { echo "Error in TOBIAS FootprintScores step for control"; exit 1; }
fi

echo -e "\nRunning TOBIAS BINDetect on treatment sample to identify bound motifs\n"

# BINDetect for treatment sample
singularity exec --bind "$PWD" \
"$TOBIAS_CONTAINER" \
TOBIAS BINDetect \
--motifs "$motifs" \
--signal "${treat_prefix}_results/${treat_prefix}_corrected.bw" \
--regions "$treat_peaks" \
--genome "$genome_fasta" \
--output "${treat_prefix}_results/bindetect_treat" || { echo "Error in TOBIAS BINDetect step for treatment"; exit 1; }

# Run BINDetect for differential analysis if control is provided
if [ ! -z "$control_bam" ] && [ ! -z "$control_peaks" ]; then
  echo -e "\nRunning TOBIAS BINDetect for differential analysis between treatment and control\n"
  
  singularity exec --bind "$PWD" \
  "$TOBIAS_CONTAINER" \
  TOBIAS BINDetect \
  --motifs "$motifs" \
  --signals "${treat_prefix}_results/${treat_prefix}_corrected.bw" "${treat_prefix}_results/${control_prefix}_corrected.bw" \
  --regions "$treat_peaks" "$control_peaks" \
  --genome "$genome_fasta" \
  --output "${treat_prefix}_results/bindetect_differential" \
  --cond_names "treatment" "control" || { echo "Error in TOBIAS BINDetect differential analysis"; exit 1; }
fi

echo -e "\nTOBIAS analysis completed successfully\n"
echo -e "Results are available in ${treat_prefix}_results directory\n" 