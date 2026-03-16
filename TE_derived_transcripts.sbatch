#!/bin/bash
#SBATCH --job-name=te_domain_overlap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# =============================================================
# Configuration
# =============================================================
MISSING_DOMAINS=/path/to/collapsed_truncations.txt
TE_BOUNDARY=/path/to/te_boundary_intersect.txt
TRANSDECODER_GFF=/path/to/transdecoder.genome.gff3
EXPRESSION=/path/to/filtered_transcripts.txt
IMMUNE_IDS=/path/to/merged_unique_ids.txt
OUTDIR=/path/to/te_domain_overlap/results

# =============================================================
# Job info
# =============================================================
echo "Job name:       $SLURM_JOB_NAME"
echo "Job ID:         $SLURM_JOB_ID"
echo "Nodes:          $SLURM_JOB_NODELIST"
echo "Date:           $(date)"

# =============================================================
# Setup
# =============================================================
mkdir -p ${OUTDIR}

# =============================================================
# Step 1: Subset missing domains to immune genes
# =============================================================
echo "Subsetting missing domains to immune genes..."
awk '
BEGIN { FS=OFS="\t" }
FNR==NR { immune[$1]; next }
($2 in immune) { print $0 }
' ${IMMUNE_IDS} ${MISSING_DOMAINS} > ${OUTDIR}/immune_truncations_full.txt

echo "Immune truncated transcripts: $(wc -l < ${OUTDIR}/immune_truncations_full.txt)"

# =============================================================
# Step 2: Build CDS exon chain per transcript from TransDecoder
# =============================================================
echo "Building CDS exon chain from TransDecoder GFF3..."
awk '
BEGIN { FS=OFS="\t" }
$3 == "CDS" {
    match($9, /Parent=([^;]+)/, a); tid=a[1]
    gsub(/\.p[0-9]*$/, "", tid)
    chr[tid] = $1
    strand[tid] = $7
    exon_count[tid]++
    exons[tid][exon_count[tid]] = $4"\t"$5
}
END {
    for(tid in exon_count) {
        n = exon_count[tid]
        for(i=1;i<n;i++) for(j=i+1;j<=n;j++) {
            split(exons[tid][i],a,"\t"); split(exons[tid][j],b,"\t")
            if(a[1]+0 > b[1]+0) { tmp=exons[tid][i]; exons[tid][i]=exons[tid][j]; exons[tid][j]=tmp }
        }
        if(strand[tid] == "-") {
            for(i=1;i<=int(n/2);i++) {
                tmp=exons[tid][i]; exons[tid][i]=exons[tid][n-i+1]; exons[tid][n-i+1]=tmp
            }
        }
        cum=0
        for(i=1;i<=n;i++) {
            split(exons[tid][i],e,"\t")
            exon_len = e[2]-e[1]+1
            print tid, chr[tid], e[1], e[2], strand[tid], cum, cum+exon_len-1
            cum += exon_len
        }
    }
}
' ${TRANSDECODER_GFF} > ${OUTDIR}/cds_exon_chain.txt

echo "CDS exon chain entries: $(wc -l < ${OUTDIR}/cds_exon_chain.txt)"
head -3 ${OUTDIR}/cds_exon_chain.txt

# =============================================================
# Step 3: Convert truncation AA coords to genomic coords
# =============================================================
echo "Converting truncation AA coordinates to genomic coordinates..."
awk '
BEGIN { FS=OFS="\t" }
FNR==NR {
    tid=$1
    chain[tid][++chain_count[tid]] = $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7
    next
}
{
    tid=$1; trunc_aa_start=$3+0; trunc_aa_end=$4+0
    trunc_nt_start = trunc_aa_start * 3
    trunc_nt_end   = trunc_aa_end * 3

    if(!(tid in chain)) next

    n = chain_count[tid]
    gstart=-1; gend=-1
    gchr=""; gstrand=""

    for(i=1;i<=n;i++) {
        split(chain[tid][i], e, "\t")
        echr=e[1]; estart=e[2]+0; eend=e[3]+0; strand=e[4]
        cum_start=e[5]+0; cum_end=e[6]+0

        if(gstart==-1 && trunc_nt_start >= cum_start && trunc_nt_start <= cum_end) {
            offset = trunc_nt_start - cum_start
            if(strand == "+") gstart = estart + offset
            else              gstart = eend - offset
            gchr = echr
            gstrand = strand
        }
        if(gend==-1 && trunc_nt_end >= cum_start && trunc_nt_end <= cum_end) {
            offset = trunc_nt_end - cum_start
            if(strand == "+") gend = estart + offset
            else              gend = eend - offset
        }
    }

    if(gstart != -1 && gend == -1) {
        split(chain[tid][n], e, "\t")
        if(gstrand == "+") gend = e[3]+0
        else               gend = e[2]+0
    }

    if(gstart != -1 && gend != -1) {
        if(gstart > gend) { tmp=gstart; gstart=gend; gend=tmp }
        print tid, gchr, gstart, gend, $0
    }
}
' ${OUTDIR}/cds_exon_chain.txt ${OUTDIR}/immune_truncations_full.txt > ${OUTDIR}/immune_truncations_genomic.txt

echo "Truncations with genomic coords: $(wc -l < ${OUTDIR}/immune_truncations_genomic.txt)"
head -3 ${OUTDIR}/immune_truncations_genomic.txt

# =============================================================
# Step 4: Internal CDS hits — TE boundary within truncation region
# =============================================================
echo "Identifying internal CDS hits..."
awk '
BEGIN { FS=OFS="\t" }
FNR==NR {
    tid=$4; gsub(/\.p[0-9]*$/, "", tid)
    if($5 != "internal_donor" && $5 != "internal_acceptor") next
    te_info=$1":"$2"-"$3";"$5";"$10";"$11
    te[tid][++te_count[tid]] = $1"\t"($2+0)"\t"te_info
    next
}
{
    tid=$1; gchr=$2; gstart=$3+0; gend=$4+0
    if(!(tid in te)) next
    for(i=1;i<=te_count[tid];i++) {
        split(te[tid][i], t, "\t")
        if(t[1]==gchr && t[2]+0 >= gstart && t[2]+0 <= gend) {
            print $0, t[3], "internal_CDS"
        }
    }
}
' ${TE_BOUNDARY} ${OUTDIR}/immune_truncations_genomic.txt > ${OUTDIR}/hits_internal_cds.txt

echo "Internal CDS hits: $(wc -l < ${OUTDIR}/hits_internal_cds.txt)"

# =============================================================
# Step 5: External UTR hits — TE within 100bp of outermost
#         boundary AND >=10% exon overlap
# =============================================================
# TE_BOUNDARY columns:
#  1=chr  2=exon_start  3=exon_end  4=tid  5=boundary_type  6=transcript_strand
#  7=TE_chr  8=TE_start  9=TE_end  10=TE_name  11=milliDiv  12=TE_strand  13=overlap_bp
echo "Identifying external UTR hits..."
awk '
BEGIN { FS=OFS="\t" }
FNR==NR {
    tid=$4; gsub(/\.p[0-9]*$/, "", tid)
    if($5 != "external_5prime" && $5 != "external_3prime") next

    exon_start=$2+0; exon_end=$3+0
    te_start=$8+0;   te_end=$9+0
    strand=$6
    btype=$5
    overlap=$13+0
    exon_len = exon_end - exon_start

    # Outermost boundary coordinate (strand-aware)
    if(btype == "external_5prime") {
        outer = (strand == "+") ? exon_start : exon_end
    } else {
        outer = (strand == "+") ? exon_end : exon_start
    }

    # Filter 1: TE within 100bp of outermost boundary
    if(te_end < outer)        te_dist = outer - te_end
    else if(te_start > outer) te_dist = te_start - outer
    else                      te_dist = 0
    if(te_dist > 100) next

    # Filter 2: TE overlaps terminal exon by >=10%
    if(exon_len == 0) next
    if((overlap / exon_len) < 0.10) next

    utr_label = (btype == "external_5prime") ? "UTR_5prime" : "UTR_3prime"
    te_info=$1":"$2"-"$3";"btype";"$10";"$11
    te[tid][++te_count[tid]] = te_info"\t"utr_label
    next
}
{
    tid=$1
    if(!(tid in te)) next
    for(i=1;i<=te_count[tid];i++) {
        split(te[tid][i], t, "\t")
        print $0, t[1], t[2]
    }
}
' ${TE_BOUNDARY} ${OUTDIR}/immune_truncations_full.txt > ${OUTDIR}/hits_external_utr.txt

echo "External UTR hits: $(wc -l < ${OUTDIR}/hits_external_utr.txt)"

# =============================================================
# Step 6: Add expression info to internal CDS hits
# =============================================================
echo "Adding expression info to internal CDS hits..."
awk '
BEGIN { FS=OFS="\t" }
FNR==NR {
    if(NR==1) next
    expr[$1] = $2"\t"$3"\t"$4"\t"$5
    next
}
{ if($1 in expr) print $0, expr[$1] }
' ${EXPRESSION} ${OUTDIR}/hits_internal_cds.txt > ${OUTDIR}/hits_internal_cds_expr.txt

echo "Internal CDS hits with expression: $(wc -l < ${OUTDIR}/hits_internal_cds_expr.txt)"

# =============================================================
# Step 7: Add expression info to external UTR hits
# =============================================================
echo "Adding expression info to external UTR hits..."
awk '
BEGIN { FS=OFS="\t" }
FNR==NR {
    if(NR==1) next
    expr[$1] = $2"\t"$3"\t"$4"\t"$5
    next
}
{ if($1 in expr) print $0, expr[$1] }
' ${EXPRESSION} ${OUTDIR}/hits_external_utr.txt > ${OUTDIR}/hits_external_utr_expr.txt

echo "External UTR hits with expression: $(wc -l < ${OUTDIR}/hits_external_utr_expr.txt)"

# =============================================================
# Step 8: Write final tables
# =============================================================
echo "Writing final tables..."

# Internal CDS final table
echo -e "transcript_id\tuniprot_id\tgene_name\ttrunc_aa_start\ttrunc_aa_end\ttrunc_genomic_coords\tgene_id\ttissue\tCPM\tisoform_fraction\tte_boundary_info\tte_location" \
    > ${OUTDIR}/final_internal_cds.txt

awk '
BEGIN { FS=OFS="\t" }
{
    tid=$1; gchr=$2; gstart=$3; gend=$4
    uniprot=$6; trunc_start=$7; trunc_end=$8
    trunc_coords=gchr":"gstart"-"gend
    gene_name=$14
    te_info=$(NF-5)
    te_location=$(NF-4)
    gene_id=$(NF-3)
    tissue=$(NF-2)
    cpm=$(NF-1)
    isoform_frac=$NF
    print tid, uniprot, gene_name, trunc_start, trunc_end, trunc_coords, gene_id, tissue, cpm, isoform_frac, te_info, te_location
}
' ${OUTDIR}/hits_internal_cds_expr.txt >> ${OUTDIR}/final_internal_cds.txt

echo "Internal CDS final table: $(wc -l < ${OUTDIR}/final_internal_cds.txt) lines (including header)"

# External UTR final table
echo -e "transcript_id\tuniprot_id\tgene_name\ttrunc_aa_start\ttrunc_aa_end\tgene_id\ttissue\tCPM\tisoform_fraction\tte_boundary_info\tte_location" \
    > ${OUTDIR}/final_external_utr.txt

awk '
BEGIN { FS=OFS="\t" }
{
    tid=$1
    uniprot=$2; trunc_start=$3; trunc_end=$4
    gene_name=$10
    te_info=$(NF-5)
    te_location=$(NF-4)
    gene_id=$(NF-3)
    tissue=$(NF-2)
    cpm=$(NF-1)
    isoform_frac=$NF
    print tid, uniprot, gene_name, trunc_start, trunc_end, gene_id, tissue, cpm, isoform_frac, te_info, te_location
}
' ${OUTDIR}/hits_external_utr_expr.txt >> ${OUTDIR}/final_external_utr.txt

echo "External UTR final table: $(wc -l < ${OUTDIR}/final_external_utr.txt) lines (including header)"

# =============================================================
# Summary
# =============================================================
echo ""
echo "=== Summary ==="
echo "Internal CDS hits with expression: $(wc -l < ${OUTDIR}/hits_internal_cds_expr.txt)"
echo "External UTR hits with expression: $(wc -l < ${OUTDIR}/hits_external_utr_expr.txt)"
echo ""
echo "TE location breakdown - Internal:"
awk 'NR>1 {print $NF}' ${OUTDIR}/final_internal_cds.txt | sort | uniq -c
echo ""
echo "TE location breakdown - External:"
awk 'NR>1 {print $NF}' ${OUTDIR}/final_external_utr.txt | sort | uniq -c
echo ""
echo "Output files:"
echo "  ${OUTDIR}/final_internal_cds.txt"
echo "  ${OUTDIR}/final_external_utr.txt"
echo "Done: $(date)"
