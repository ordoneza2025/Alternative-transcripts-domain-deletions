#!/bin/bash
#SBATCH --job-name=transdecoder_predict
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
BLASTP_OUT=/path/to/blastp.outfmt6
OUTDIR=/path/to/transdecoder/predict
GENOME_GFF=/path/to/genome.gff3        # genome-based annotation (for util step)
GENOME_FA=/path/to/genome.fa           # genome FASTA (for util step)

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
echo "TransDecoder version:"; TransDecoder.Predict --version

# =============================================================
# Setup
# =============================================================
mkdir -p ${OUTDIR}
cd ${OUTDIR}

# =============================================================
# Step 1: Predict coding regions (with blastp support)
# =============================================================
echo "Running TransDecoder.Predict with blastp evidence..."
TransDecoder.Predict \
    -t ${TRANSCRIPTS} \
    --retain_blastp_hits ${BLASTP_OUT} \
    --single_best_only

# =============================================================
# Step 2: Generate genome-based coordinates
# =============================================================
echo "Generating genome-based GFF3..."
TransDecoder.Util.gtf_to_alignment_gff3.pl ${GENOME_GFF} > genome_aligned.gff3

cdna_alignment_orf_to_genome_orf.pl \
    ${OUTDIR}/$(basename ${TRANSCRIPTS}).transdecoder.gff3 \
    genome_aligned.gff3 \
    ${TRANSCRIPTS} \
    > $(basename ${TRANSCRIPTS}).transdecoder.genome.gff3

# =============================================================
# Step 3: Summary of output files
# =============================================================
PREFIX=${OUTDIR}/$(basename ${TRANSCRIPTS}).transdecoder

echo "Output files:"
for ext in pep cds gff3 bed; do
    f="${PREFIX}.${ext}"
    if [[ -f "$f" ]]; then
        echo "  $f ($(wc -l < $f) lines)"
    fi
done

echo "Predicted ORF count:"
grep -c ">" ${PREFIX}.pep

echo "Done: $(date)"
