#!/bin/bash
# Submit the 26 assembly tasks in parallel, then the summary job, which waits for all of them.
# Run from /N/scratch/vgutta/assembly:   bash submit_all.sh
cd "$(dirname "$0")"
mkdir -p logs
ARRAY_ID=$(sbatch --parsable assemble_array.sh)
SUM_ID=$(sbatch --parsable --dependency=afterany:${ARRAY_ID} summarize.sh)
echo "Assembly array job: $ARRAY_ID  (26 tasks)"
echo "Summary job:        $SUM_ID  (starts when all tasks finish)"
echo "Check progress:     squeue -u $USER"
