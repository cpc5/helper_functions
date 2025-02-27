#!/bin/bash

# Parse command-line options
while getopts n:h arg; do
  case "$arg" in
    n) n="$OPTARG" ;;
    h)
      echo "Usage: FE_bdg_to_bigwig.sh -n <sample_prefix>"
      echo "   -n: Sample name prefix (e.g., mESC_SUZ12 for mESC_SUZ12_control_lambda.bdg)"
      echo "   -h: Display this help message"
      exit 0
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# Check if sample name is provided
if [ -z "$n" ]; then
  echo "Error: sample name prefix not provided. Use -n to specify the sample name prefix"
  exit 1
fi

module load singularity

# Cleanup temp files on exit
trap "rm -f xtmp xtmp_FE.bdg.clip" EXIT

echo "Running bdg compare\n"

singularity exec --bind "$PWD" \
/project/imoskowitz/shared/software/singularityImages/singularity-images/depot.galaxyproject.org-singularity-macs2-2.2.7.1--py38h4a8c8d9_3.img \
macs2 bdgcmp -t "${n}_treat_pileup.bdg" \
        -c "${n}_control_lambda.bdg" \
        -o "${n}_FE.bdg" -m FE || { echo "Error in MACS2 step"; exit 1; }

echo "Running bedtools slop\n"

singularity exec --bind "$PWD" \
/project/imoskowitz/shared/software/singularityImages/singularity-images/depot.galaxyproject.org-singularity-bedtools-2.30.0--hc088bd4_0.img \
bedtools slop -i \
"${n}_FE.bdg" \
-g /project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Sequence/mm10.chrom.sizes \
-b 0 > xtmp || { echo "Error in bedtools slop step"; exit 1; }

echo "Running bedClip to Remove lines from bed file that refer to off-chromosome locations\n"

singularity exec --bind "$PWD" \
/project/imoskowitz/shared/software/singularityImages/singularity-images/depot.galaxyproject.org-singularity-ucsc-bedclip-377--h0b8a92a_2.img \
bedClip xtmp /project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Sequence/mm10.chrom.sizes \
xtmp_FE.bdg.clip || { echo "Error in bedClip step"; exit 1; }


echo "sorting FE bdg.."
LC_COLLATE=C sort -k1,1 -k2,2n xtmp_FE.bdg.clip > xtmp

echo "converting bdg to bw "
singularity exec --bind "$PWD" \
/project/imoskowitz/shared/software/singularityImages/singularity-images/depot.galaxyproject.org-singularity-ucsc-bedgraphtobigwig-377--h446ed27_1.img \
bedGraphToBigWig \
xtmp \
/project/imoskowitz/shared/annotations/iGenomes/Mus_musculus/UCSC/mm10/Sequence/mm10.chrom.sizes \
"${n}_FE.bigWig" || { echo "Error in bedGraphToBigWig step"; exit 1; }

# Temporary files will be cleaned up on exit

