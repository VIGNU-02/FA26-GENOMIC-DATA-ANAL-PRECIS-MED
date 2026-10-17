#!/bin/bash
#SBATCH -J ecoli_asm
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-25
#SBATCH -o logs/asm_%A_%a.log
#
# Job array: one task per (tool, k). Tasks 0-12 = Velvet, 13-25 = Oases, same 13 k values.
# If you change K, also change --array to 0-(2*number_of_k - 1).
# Submit with submit_all.sh (it also queues the summary job).

SRR="SRR21904868"
CONDA_ENV="info636"
K=(21 31 41 51 61 71 81 91 101 111 121 131 141)

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
mkdir -p "$WORK"/{velvet,oases} "$ROWS"

# ---------- which tool and k is this task? ----------
i=${SLURM_ARRAY_TASK_ID:?run this with sbatch}
n=${#K[@]}
if (( i < n )); then tool=velvet; else tool=oases; fi
k=${K[$(( i % n ))]}
d="$WORK/$tool/k$k"
echo "Task $i: $tool k=$k on $(hostname), $THREADS threads"

# ---------- reads (already downloaded with prefetch + fasterq-dump --split-3) ----------
if [[ -f $READS/${SRR}_1.fastq && -f $READS/${SRR}_2.fastq ]]; then
    READ_ARGS=(-shortPaired -fastq -separate "$READS/${SRR}_1.fastq" "$READS/${SRR}_2.fastq")
else
    echo "Reads not found in $READS" >&2; exit 1
fi

# time each command: wall seconds + peak memory (KB); wall time only if GNU time is missing
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

# ---------- assemble ----------
rm -rf "$d"
if [[ $tool == velvet ]]; then
    timed "$d.velveth.time" velveth "$d" "$k" "${READ_ARGS[@]}"
    timed "$d.velvetg.time" velvetg "$d" -exp_cov auto -cov_cutoff auto
    fasta="$d/contigs.fa"
    times=("$d.velveth.time" "$d.velvetg.time")
else
    timed "$d.velveth.time" velveth "$d" "$k" "${READ_ARGS[@]}"
    timed "$d.velvetg.time" velvetg "$d" -read_trkg yes
    timed "$d.oases.time"   oases "$d"
    fasta="$d/transcripts.fa"
    times=("$d.velveth.time" "$d.velvetg.time" "$d.oases.time")
fi

# ---------- one summary row for this assembly ----------
read wall mem < <(combine "${times[@]}")
python3 "$STATS" summary "$tool" "$k" "$fasta" "$wall" "$mem" > "$ROWS/${tool}_k${k}.tsv"
cat "$ROWS/${tool}_k${k}.tsv"

# Velvet's big intermediate files are no longer needed once contigs exist
rm -f "$d/Sequences" "$d/Roadmaps" "$d/PreGraph"
echo "Task $i done"
