#!/bin/bash
# submit a job to SLURM quickly. by Carlos Perez-Cervantes
# Initialize variables with default values
n="myrun"
c="8"
f="echo Hello world"
m="32"
t="04"
l=""
p="caslake"

# Parse command-line options
while getopts n:c:p:f:m:l:t:h arg; do
  case "$arg" in
    n) n="$OPTARG" ;;
    c) c="$OPTARG" ;;
    f) f="$OPTARG" ;;
    m) m="$OPTARG" ;;
    t) t="$OPTARG" ;;
    l) l="$OPTARG" ;;
    p) p="$OPTARG" ;;
    h)
      echo "Usage:"
      echo "   quickslurm.sh -h   Display this help message."
      echo "   -n: Job name"
      echo "   -c: Number of CPU cores per task"
      echo "   -f: Specify command to execute"
      echo "   -m: Memory per node in GB"
      echo "   -t: Walltime in hours"
      echo "   -l: Module"
      echo "   -p: Partition/queue name"
      exit 0
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# Check if required options are provided
#if ! [ -z "$n" ] || [ -z "$c" ] || [ -z "$f" ] || [ -z "$m" ] || [ -z "$t" ] || [ -z "$l" ] || [ -z "$p" ]; then
#  echo "Missing required options. Use -h for help."
#  exit 1
#fi

# Write SLURM script
cat << EOM > "$n.sh"
#!/bin/bash
#SBATCH --job-name=$n
#SBATCH --output=$n.log
#SBATCH --error=$n.err
#SBATCH --partition=$p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=$c
#SBATCH --mem=${m}GB
#SBATCH --time=$t:00:00
#SBATCH --account=pi-imoskowitz

$l

cd $(pwd)

$f

EOM

echo "SLURM script '$n.sh' created."
