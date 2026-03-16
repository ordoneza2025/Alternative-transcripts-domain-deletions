#!/bin/bash
#SBATCH --job-name=te_cds_boundary
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# =============================================================
# Configuration
# =============================================================
CDS_GFF=/path/to/transdecoder_cds.gff3
TRANSCRIPT_GTF=/path/to/transcripts.gtf
REPEATS=/path/to/repeats.rmsk.bed
OUTDIR=/path/to/te_boundary_results_CDS

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
module load bedtools
echo "bedtools version:"; bedtools --version | head -1

# =============================================================
# Setup
# =============================================================
mkdir -p ${OUTDIR}

echo "=== Inputs ==="
echo "CDS GFF3:       ${CDS_GFF}"
echo "Transcript GTF: ${TRANSCRIPT_GTF}"
echo "Repeats:        ${REPEATS}"
echo "Output:         ${OUTDIR}"

# =============================================================
# Step 1: Get full exon order from transcript GTF
# =============================================================
echo "Getting full exon order from transcript GTF..."
awk '
BEGIN { FS=OFS="\t" }
$3 == "exon" {
    tid=""
    if(match($0, /transcript_id "([^"]+)"/, a)) tid=a[1]
    if(tid == "") next
    exon_count[tid]++
    exons[tid][exon_count[tid]] = $4"\t"$5"\t"$7
}
END {
    for(tid in exon_count) {
        n = exon_count[tid]
        for(i=1;i<n;i++) for(j=i+1;j<=n;j++) {
            split(exons[tid][i],a,"\t"); split(exons[tid][j],b,"\t")
            if(a[1]+0 > b[1]+0) { tmp=exons[tid][i]; exons[tid][i]=exons[tid][j]; exons[tid][j]=tmp }
        }
        for(i=1;i<=n;i++) {
            split(exons[tid][i], e, "\t")
            print tid, e[1], e[2], e[3], i, n
        }
    }
}
' ${TRANSCRIPT_GTF} > ${OUTDIR}/transcript_exon_order.txt

echo "Transcript exon records: $(wc -l < ${OUTDIR}/transcript_exon_order.txt)"
head -3 ${OUTDIR}/transcript_exon_order.txt

# =============================================================
# Step 2: Extract CDS boundary points using transcript exon
#         order for labeling
# =============================================================
echo "Extracting CDS boundary points..."
awk '
BEGIN { FS=OFS="\t" }
FNR==NR {
    key=$1"\t"$2"\t"$3
    rank[$1][key] = $5
    total[$1] = $6
    next
}
$3 == "CDS" {
    chr=$1; start=$4-1; end=$5; strand=$7
    tid=""
    if(match($9, /Parent=([^;]+)/, a)) tid=a[1]
    if(tid == "") next
    base_tid = tid; gsub(/\.p[0-9]*$/, "", base_tid)

    found_rank = -1; found_total = -1
    if(base_tid in rank) {
        found_total = total[base_tid]
        for(key in rank[base_tid]) {
            split(key, k, "\t")
            flair_start = k[2]+0; flair_end = k[3]+0
            if(start+0 >= flair_start && end+0 <= flair_end+1) {
                found_rank = rank[base_tid][key]
                break
            }
        }
    }

    # Account for minus strand rank reversal
    if(strand == "-") {
        is_first = (found_rank == found_total)
        is_last  = (found_rank == 1)
    } else {
        is_first = (found_rank == 1)
        is_last  = (found_rank == found_total)
    }
    if(found_rank == -1) { is_first=0; is_last=0 }

    if(strand == "+") {
        if(is_first) print chr, start,   start+1, tid, "external_5prime",   strand
        else          print chr, start,   start+1, tid, "internal_acceptor", strand
        if(is_last)  print chr, end-1,   end,     tid, "external_3prime",   strand
        else          print chr, end-1,   end,     tid, "internal_donor",    strand
    } else {
        if(is_first) print chr, end-1,   end,     tid, "external_5prime",   strand
        else          print chr, end-1,   end,     tid, "internal_acceptor", strand
        if(is_last)  print chr, start,   start+1, tid, "external_3prime",   strand
        else          print chr, start,   start+1, tid, "internal_donor",    strand
    }
}
' ${OUTDIR}/transcript_exon_order.txt ${CDS_GFF} > ${OUTDIR}/cds_boundaries.bed

echo "Boundary type counts:"
cut -f5 ${OUTDIR}/cds_boundaries.bed | sort | uniq -c
echo "Sample:"; head -3 ${OUTDIR}/cds_boundaries.bed

# =============================================================
# Step 3: Intersect CDS boundaries with TEs
# =============================================================
echo "Intersecting boundaries with TEs..."
bedtools intersect \
    -a ${OUTDIR}/cds_boundaries.bed \
    -b ${REPEATS} \
    -wo \
    > ${OUTDIR}/te_boundary_intersect.txt

echo "TE-boundary intersections: $(wc -l < ${OUTDIR}/te_boundary_intersect.txt)"
echo "Columns in output: $(head -1 ${OUTDIR}/te_boundary_intersect.txt | awk '{print NF}')"
echo "Sample:"; head -3 ${OUTDIR}/te_boundary_intersect.txt

# =============================================================
# Step 4: Summarize by boundary type
# =============================================================
echo "Summarizing by boundary type..."
awk 'BEGIN{OFS="\t"} {print $4, $5, $10}' ${OUTDIR}/te_boundary_intersect.txt \
    | sort | uniq \
    > ${OUTDIR}/te_boundary_summary.txt

echo "Boundary type counts:"
cut -f2 ${OUTDIR}/te_boundary_summary.txt | sort | uniq -c

echo ""
echo "Done. Output files:"
echo "  ${OUTDIR}/cds_boundaries.bed"
echo "  ${OUTDIR}/te_boundary_intersect.txt"
echo "  ${OUTDIR}/te_boundary_summary.txt"
echo "Done: $(date)"
