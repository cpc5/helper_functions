#!/bin/bash

### Merge BAM files using samtools ###
### Carlos Perez-Cervantes Moskowitz lab

# Parse command-line options
while getopts o:h arg; do
  case "$arg" in
    o) output="$OPTARG" ;;
    h)
      echo "Usage: merge_bams.sh -o <output_bam> input1.bam input2.bam [input3.bam ...]"
      echo "   -o: Output BAM file name"
      echo "   -h: Display this help message"
      echo ""
      echo "All input BAM files should be listed after the options"
      exit 0
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# Shift past the options to get to the input files
shift $((OPTIND-1))

# Check if output name is provided
if [ -z "$output" ]; then
  echo "Error: output BAM file name not provided"
  echo "Use -o to specify the output BAM file name"
  exit 1
fi

# Check if at least two input files are provided
if [ $# -lt 2 ]; then
  echo "Error: at least two input BAM files are required"
  echo "Usage: merge_bams.sh -o <output_bam> input1.bam input2.bam [input3.bam ...]"
  exit 1
fi

# Check if all input files exist
for file in "$@"; do
  if [ ! -f "$file" ]; then
    echo "Error: input file '$file' does not exist"
    exit 1
  fi
done

module load singularity

# Define the samtools container
SAMTOOLS_CONTAINER="/project/imoskowitz/shared/software/singularityImages/singularity-images/depot.galaxyproject.org-singularity-samtools-1.17--h00cdaf9_0.img"

echo -e "\nMerging BAM files...\n"
echo "Input files:"
for file in "$@"; do
  echo "  - $file"
done
echo "Output: $output"

# Run samtools merge
singularity exec --bind "$PWD" \
"$SAMTOOLS_CONTAINER" \
samtools merge -o "$output" "$@" || { echo "Error in samtools merge step"; exit 1; }

echo -e "\nBAM merge completed successfully\n"
echo "Output BAM file: $output" 