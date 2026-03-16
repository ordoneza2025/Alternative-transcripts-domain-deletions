#!/bin/bash
#SBATCH --job-name=transdecoder_longorfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --partition=short
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# =============================================================
# Configuration
# =============================================================
TRANSCRIPTS=/path/to/transcripts.fa
OUTDIR=/path/to/transdecoder/longorfs
MIN_AA=50

# =============================================================
# Job info
# =============================================================
echo "Job name:       $SLURM_JOB_NAME"
echo "Job ID:         $SLURM_JOB_ID"
echo "Nodes:          $SLURM_JOB_NODELIST"
echo "Date:           $(date)"

# =============================================================
# Load modules
# =============================================================
module load transdecoder
echo "TransDecoder version:"; TransDecoder.LongOrfs --version

# =============================================================
# Setup
# =============================================================
mkdir -p ${OUTDIR}
cd ${OUTDIR}

# =============================================================
# Step 1: Identify long ORFs
# =============================================================
echo "Identifying long ORFs (min AA length: ${MIN_AA})..."
TransDecoder.LongOrfs \
    -t ${TRANSCRIPTS} \
    -m ${MIN_AA}

# =============================================================
# Step 2: Summary
# =============================================================
echo "Output directory contents:"
ls -lh ${OUTDIR}/$(basename ${TRANSCRIPTS}).transdecoder_dir/

echo "Candidate ORF count:"
grep -c ">" ${OUTDIR}/$(basename ${TRANSCRIPTS}).transdecoder_dir/longest_orfs.pep

echo "Done: $(date)"
