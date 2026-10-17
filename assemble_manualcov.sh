#!/bin/bash
#SBATCH -J ecoli_mcov
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-5
#SBATCH -o logs/mcov_%A_%a.log
# Velvet at small k with expected coverage set manually instead of "-exp_cov auto".
# At k <= 71 velvet's auto estimate failed (e.g. k=51 estimated coverage 1.0), so no
# error k-mers were removed. Here exp_cov is scaled from the k=91 estimate (26.6):
#   exp_cov(k) = 26.6 * (L - k + 1) / (L - 91 + 1),  L = 151,   cov_cutoff = exp_cov / 2

SRR="SRR21904868"
CONDA_ENV="info636"
K=(21 31 41 51 61 71)

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
STATS="$SUBMIT_DIR/assembly_stats.py"
WORK="$SUBMIT_DIR/ecoli_assembly"
READS="$WORK/reads"
ROWS="$WORK/results/rows"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="$THREADS" OMP_THREAD_LIMIT="$THREADS"

module load conda
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"
set -eo pipefail
mkdir -p "$WORK/velvet_manualcov" "$ROWS"

k=${K[${SLURM_ARRAY_TASK_ID:?run this with sbatch}]}
EXP=$(awk -v k="$k" 'BEGIN{printf "%.1f", 26.6*(152-k)/61}')
CUT=$(awk -v e="$EXP" 'BEGIN{printf "%.1f", e/2}')
d="$WORK/velvet_manualcov/k$k"
echo "velvet (manual coverage) k=$k exp_cov=$EXP cov_cutoff=$CUT on $(hostname)"

timed() {
    local out=$1; shift
    if [[ -x /usr/bin/time ]]; then
        /usr/bin/time -f "%e\t%M" -o "$out" "$@"
    else
        local t0=$(date +%s.%N); "$@"; local t1=$(date +%s.%N)
        awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f\t0\n", b-a}' > "$out"
    fi
}
combine() { awk -F'\t' 'NF==2{w+=$1; if($2>m)m=$2} END{printf "%.1f %d\n", w, m}' "$@"; }

rm -rf "$d"
timed "$d.velveth.time" velveth "$d" "$k" -shortPaired -fastq -separate \
    "$READS/${SRR}_1.fastq" "$READS/${SRR}_2.fastq"
timed "$d.velvetg.time" velvetg "$d" -exp_cov "$EXP" -cov_cutoff "$CUT"

read wall mem < <(combine "$d.velveth.time" "$d.velvetg.time")
python3 "$STATS" summary velvet_manualcov "$k" "$d/contigs.fa" "$wall" "$mem" > "$ROWS/velvetmcov_k${k}.tsv"
echo -e "$k\t$EXP\t$CUT" > "$ROWS/velvetmcov_k${k}.params"
cat "$ROWS/velvetmcov_k${k}.tsv"
rm -f "$d/Sequences" "$d/Roadmaps" "$d/PreGraph"
echo done
