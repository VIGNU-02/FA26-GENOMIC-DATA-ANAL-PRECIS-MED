#!/bin/bash
#SBATCH -J ecoli_summary
#SBATCH -p general
#SBATCH -A c02374
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH -o logs/summary_%j.log
# Runs after every array task has finished: builds the summary table, picks the best k
# for each tool (highest N50) and lists every contig length for those assemblies.

SRR="SRR21904868"
CONDA_ENV="info636"

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
STATS="$SUBMIT_DIR/assembly_stats.py"
WORK="$SUBMIT_DIR/ecoli_assembly"
READS="$WORK/reads"
RES="$WORK/results"

module load conda
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"
set -eo pipefail

# read statistics
read NREADS TOTBP MAXLEN MEANLEN < <(cat "$READS/${SRR}_1.fastq" "$READS/${SRR}_2.fastq" |
    awk 'NR%4==2{n++; b+=length($0); if(length($0)>m)m=length($0)} END{printf "%d %d %d %.1f\n", n, b, m, b/n}')
{
    echo -e "accession\t$SRR"
    echo -e "layout\tpaired-end"
    echo -e "total_reads\t$NREADS"
    echo -e "total_bases\t$TOTBP"
    echo -e "mean_read_length\t$MEANLEN"
    echo -e "max_read_length\t$MAXLEN"
    echo -e "est_coverage_vs_4.6Mb\t$(awk -v b="$TOTBP" 'BEGIN{printf "%.1fx", b/4.6e6}')"
} | tee "$RES/read_info.tsv"

# summary table: velvet rows then oases rows, each sorted by k
{
    python3 "$STATS" header
    for tool in velvet oases; do
        cat "$RES"/rows/${tool}_k*.tsv 2>/dev/null | sort -t$'\t' -k2,2n
    done
} > "$RES/summary.tsv"

# Velvet's own Log files (N50 and totals as velvet reports them)
for f in "$WORK"/velvet/k*/Log; do
    k=$(basename "$(dirname "$f")"); cp "$f" "$RES/velvet_${k}_Log.txt"
done

# best k per tool and contig lengths
: > "$RES/best_k.tsv"
for tool in velvet oases; do
    best=$(awk -F'\t' -v t="$tool" '$1==t' "$RES/summary.tsv" | sort -t$'\t' -k6,6nr | head -1 | cut -f2)
    [[ -z $best ]] && { echo "no $tool results"; continue; }
    if [[ $tool == velvet ]]; then fa="$WORK/velvet/k$best/contigs.fa"; else fa="$WORK/oases/k$best/transcripts.fa"; fi
    echo -e "$tool\tbest_k_by_N50\t$best" | tee -a "$RES/best_k.tsv"
    python3 "$STATS" lengths "$fa" > "$RES/${tool}_k${best}_lengths.tsv"
    cp "$fa" "$RES/${tool}_k${best}_$(basename "$fa")"
done

echo
column -t -s$'\t' "$RES/summary.tsv" 2>/dev/null || cat "$RES/summary.tsv"
echo "Done. Results are in $RES"
