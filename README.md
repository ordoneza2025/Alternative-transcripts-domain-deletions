# Identifying domain-altering alternative isoforms across mammalian long-read transcriptomes

The purpose of this workflow is to identify alternative protein-coding isoforms with predicted domain loss or truncation using long-read transcriptomes from human, mouse, and Jamaican fruit bat (*Artibeus jamaicensis*). Starting from long-read transcriptome GTF files, this pipeline predicts coding sequences, matches species-specific peptides to human reference proteins, and identifies deletions or truncations in the encoded products. Because many alternative isoforms result in partial loss of protein sequence, the main goal is to determine which protein domains are missing or disrupted.

This workflow also evaluates whether these alternative transcript structures are associated with transposable element (TE) exonization. By intersecting transcript features with repeat annotations, the pipeline identifies candidate TE-derived isoforms that may contribute to altered protein structure and function.  

## Predicting ORFs from transcripts

Long-read GTF files for human, mouse, and Jamaican fruit bat are converted into transcript FASTA format using `gffread`. Open reading frames (ORFs) are then identified with `TransDecoder.LongOrfs`, and candidate coding sequences are annotated against the human SwissProt reference using `BLASTp`. `TransDecoder.Predict` is then used together with the BLASTp results to identify the single best coding ORF for each transcript.

## Downloading and re-formatting Uniprot domain files 

These data files are found in the downloadable data for humans (hg38) in UCSC (October 2024). 
https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/

1. **UnipDomain.bb** - annotated domains
2. **UnipLocCytopl.bb** - Cytoplasmic domains.
3. **UnipLocExtra.bb** - Extracellular domains.
4. **UnipLocSignal.bb** - signal peptide sequences. 
5. **UnipLocTransMemb.bb** - transmembrane domains.
6. **UnipInterest.bb** - Other domains of interest.

These files contain both Swissprot (validated) and Trembl (computational annotation). Only Swissprot records were used, Trembl was filtered out once files were converted and reformated (grep -v "TrEMBL" input.bed > input.final.bed).  

https://github.com/ordoneza2025/missing_domains/blob/main/fromat_Uniprot_domain_files.sh
   

### 2. Finding truncations in bat transcripts

Predicted species-specific peptide sequences are matched to their corresponding human peptide sequences using the BLASTp output. Peptide pairs are then aligned with `MUSCLE` to identify regions missing from the species-specific isoform relative to the human reference. Alignment gaps greater than 10 amino acids are classified as deletions, and nearby gaps within 15 amino acids of each other are merged to simplify downstream domain assignment. The resulting deletion coordinates are reformatted into BED files containing the peptide gap coordinates and the corresponding species-specific transcript IDs

https://github.com/ordoneza2025/missing_domains/blob/main/identifying_truncations.sbatch

## Identifying truncated or deleted domains
`bedtools` is used to intersect peptide deletion coordinates with UniProt domain annotations to identify domains that are truncated or deleted in alternative isoforms. Domain categories include extracellular, cytoplasmic, signal peptide, transmembrane, and other features of interest. To focus on likely disruptive events, intersections are filtered by minimum deletion size: greater than 50 amino acids for cytoplasmic, extracellular, and other domains of interest; greater than 20 amino acids for signal peptides; and greater than 15 amino acids for transmembrane domains. Filtered domain deletion calls are merged into a final summary table.

https://github.com/ordoneza2025/Alternative-transcripts-domain-deletions/blob/main/subset.sbatch


## TE-exonization analysis

To test whether alternative transcript structures are associated with repeat-derived sequence, repeat annotations are collected for each species. Human and mouse repeat annotations are obtained from Dfam/UCSC. For Jamaican fruit bat, UCSC repeat annotations are merged with a custom bat-specific repeat library. Repeat file can be found in Zenodo: 10.5281/zenodo.19022336

Exons from significantly expressed transcripts are classified as internal or external, and exon boundary coordinates are extracted. Internal exon boundaries are intersected with repeat annotations using a 1-bp overlap criterion. External exon boundaries are intersected using a 100-bp boundary window and minimum overlap threshold of 10% to detect repeat-associated transcript start or end features.

The GFF3 input from transdecoder helps annotate the exons are internal or external. 

https://github.com/ordoneza2025/Alternative-transcripts-domain-deletions/blob/main/te_cds_boundary.sbatch

## Linking TE-derived transcript features to domain loss

To identify transcripts in which TE-associated transcript features correspond to deleted protein domains, `TransDecoder` coding annotations are used to project peptide deletion coordinates back to genomic coordinates. These projected intervals are then compared with TE-intersection coordinates to identify candidate TE-derived alternative transcripts with predicted coding consequences. Additional filters of expression and gene ontology (immune genes) are used to subset for final TE-derived transcript candidates. 

https://github.com/ordoneza2025/Alternative-transcripts-domain-deletions/blob/main/TE_derived_transcripts.sbatch

## Output

The final outputs of this workflow include:

- predicted coding ORFs for long-read transcripts
- peptide alignments between species-specific isoforms and human reference proteins
- merged deletion coordinates representing missing peptide regions
- domain intersection tables identifying truncated or deleted protein domains
- TE-exonization calls linked to alternative transcript structures
- a final set of candidate TE-derived alternative isoforms with predicted domain-altering effects







