# BRAF_analysis
Analysing BRAF V600 mutations in neoadjuvant RNA

### BRAF V600 mutation tracking

After variant calling, the filtered VCFs are used for longitudinal BRAF V600 mutation analysis (`braf_var_analysis.ipynb`). This notebook tracks BRAF V600 status and variant allele frequency (VAF) across treatment timepoints (Baseline → Week 1 → Week 2 → Week 6 → Progressed).

### V600 detection logic

BRAF V600 variants are identified at hg38 coordinates chr7:140,753,336–140,753,337 from `PASS`-filtered VCF records. The detection handles single-nucleotide and multi-nucleotide variant representations:

| Coordinate(s) | REF → ALT | Call |
|----------------|-----------|------|
| 140753336 | A → T | V600E |
| 140753336–337 | AC → TT | V600K |
| 140753336–337 | AC → CT | V600R |
| No PASS variant at locus | — | not detected |

VAF is calculated from the AD (allelic depth) field in the VCF genotype column.

### Analyses performed

1. **BRAF V600 heatmap** per-patient status (WT(not detected) / V600E / V600K/R) across timepoints, sorted by number of V600-positive timepoints, with DNA VAF and tumour cellularity side annotations
2. **Detection rate by response** proportion of samples with a detectable V600 mutation at each timepoint, stratified by MPR vs NMPR (with Wilson confidence intervals)
3. **Detection rate by treatment arm** same as above, faceted by treatment arm (e.g. PD1 vs PD1+LAG3)
4. **RNA VAF trajectories** per-patient RNA VAF over time, with group medians, stratified by response and treatment arm
5. **RNA VAF vs DNA VAF comparison** correlation of RNA-derived VAF at baseline with DNA-based VAF from Sequenza
6. **Check RNA expression** of key markers like BRAF, SOX10, MLANA as purity reference proxy

## RNA Variant Calling Pipeline

A PBS-based pipeline for RNA-seq variant calling on the NCI Gadi HPC system, following [GATK Best Practices for RNA-seq short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels).

**Author:** Jayden Beckwith

Adapted from Patrick Terrematte's pipeline for HPC of NPAD/UFRN ([original](http://hungria.imd.ufrn.br/~terrematte/aDNA/rna_seq_variant_pipeline.sh)), based on [Anand M.'s GATK RNA-seq variant pipeline](https://gist.github.com/PoisonAlien/c6c03539cf4b1ac41cf1).

---

## Overview

The pipeline takes BAMs (back to FASTQ) or FASTQ, performs STAR two-pass alignment against GRCh38, and runs the full GATK variant calling workflow. Downstream analysis includes longitudinal BRAF V600 mutation tracking across treatment timepoints. It is designed for large cohorts with built-in checkpoint/resume support and PBS job scheduling.

### Pipeline stages

1. **STAR genome index build** single job generating a STAR index from GRCh38 primary assembly + GENCODE v46 GTF
2. **FASTQ extraction + STAR two-pass alignment** per-sample: name-sort → SamToFastq → STAR `--twopassMode Basic`
3. **Variant calling** per-sample: AddOrReplaceReadGroups → MarkDuplicates → SplitNCigarReads → BaseRecalibrator → ApplyBQSR → HaplotypeCaller → VariantFiltration → VariantsToTable

All per-sample jobs depend on the index build completing successfully. Intermediate files are cleaned up automatically once all outputs are verified.

---

## Requirements

### Software (via Singularity containers)

| Tool | Container | Version |
|------|-----------|---------|
| STAR | `star_2.7.11b--h43eeafb_0.sif` | 2.7.11b |
| Picard | `picard.sif` | 3.4.0-12 |
| GATK | `gatk_4.5.0.0.sif` | 4.5.0.0 |
| Samtools | `samtools_v1.3.1_cv4.sif` | 1.3.1 |

Singularity must be available via `module load singularity` on Gadi.

### HPC

- NCI Gadi with PBS Pro scheduler
- Project allocation (default: `jo11`)
- Sufficient scratch storage at `/scratch/jo11/`

---

## Reference genome and annotation setup

The pipeline uses GRCh38 primary assembly with GENCODE v46 annotations. Below are instructions for downloading and preparing the required reference files.

### Option A: GENCODE / Ensembl primary assembly (used in this pipeline)

```bash
# Download GRCh38 primary assembly from GENCODE
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Download GENCODE v46 annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz

# Create FASTA index and sequence dictionary (needed by GATK)
samtools faidx GRCh38.primary_assembly.genome.fa

picard CreateSequenceDictionary \
    R=GRCh38.primary_assembly.genome.fa \
    O=GRCh38.primary_assembly.genome.dict
```

### Option B: Broad Institute hg38 bundle

```bash
# Download from Google Cloud genomics public data
# https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta --no-check-certificate
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai --no-check-certificate
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict --no-check-certificate
```


### Known sites for BQSR

Download these VCFs and their `.tbi` index files from the [Broad hg38 resource bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0):

- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
- `1000G_phase1.snps.high_confidence.hg38.vcf.gz`
- `dbsnp_138.hg38.vcf.gz`

> **Note:** Ensure the chromosome naming convention in your reference FASTA (e.g. `chr1` vs `1`) matches the known-sites VCFs. The GENCODE primary assembly uses `chr`-prefixed contigs, which is compatible with the Broad bundle VCFs.

---

## Directory structure

```
/scratch/jo11/
├── GRCh38.fa/
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── GRCh38.primary_assembly.genome.fa.fai
│   ├── GRCh38.primary_assembly.genome.dict
│   └── gencode.v46.chr_patch_hapl_scaff.annotation.gtf
├── STAR_index_GRCh38_gencode46/       # Built by Stage 1
├── singularity_images/
│   ├── star_2.7.11b--h43eeafb_0.sif
│   ├── picard.sif
│   ├── gatk_4.5.0.0.sif
│   └── samtools_v1.3.1_cv4.sif
├── NeoTrio_RNA/
│   ├── CAGRF220510736_rnaseq/         # Input BAMs
│   ├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
│   ├── 1000G_phase1.snps.high_confidence.hg38.vcf.gz
│   └── dbsnp_138.hg38.vcf.gz
└── NeoTrio_RNA_variants/              # Pipeline outputs
    ├── fastqs/                        # Extracted FASTQs
    ├── star_aligned/                  # STAR two-pass BAMs
    ├── pbs_logs/                      # PBS stdout/stderr
    ├── tmp/                           # Per-sample temp dirs
    └── <SAMPLE>.processed/            # Per-sample final outputs
        ├── <SAMPLE>.recal.bam         # Recalibrated BAM
        ├── <SAMPLE>.raw.vcf           # Raw HaplotypeCaller VCF
        ├── <SAMPLE>.filtered.vcf      # Filtered VCF
        └── <SAMPLE>.filtered.tsv      # Tabular variant summary
```

---

## Usage

### 1. Configure paths

Edit the `CONFIGURATION` block at the top of `RNA_variant_pipeline.sh` to set your reference paths, BAM input folder, Singularity image paths, known-sites VCFs, output directories, and PBS project/queue settings.

### 2. Submit

```bash
qsub RNA_variant_pipeline.sh
```

This submits the STAR index build job first, then submits all per-sample jobs with a dependency on the index job. Up to 30 samples run concurrently by default (configurable via `MAX_CONCURRENT`).

### 3. Skip specific samples

Set the `SKIP_SAMPLES` variable to a space-separated list of BAM basenames (without `.bam`) to exclude:

```bash
SKIP_SAMPLES="sample_A sample_B"
```

### 4. Resume after failure

Simply resubmit the pipeline script. Each step checks for its expected output before running — completed steps are skipped automatically. This means timed-out or partially failed jobs will resume from the last successful checkpoint.

---

## Pipeline details

### Checkpoint/resume

Every step uses a `checkpoint_run` function that checks whether the expected output file already exists and is non-empty before executing. If a job is killed (e.g. walltime exceeded), resubmitting the same script will resume from the last completed step.

### Submission order

Samples are sorted by timepoint before submission: PRE → ED1 → ED2 → CLND. This ensures earlier timepoints are processed first when queue capacity is limited.

### Resource allocation

| Job type | CPUs | Memory | Walltime |
|----------|------|--------|----------|
| STAR index build | 24 | 128 GB | 48 h |
| Per-sample (align + variant call) | 12 | 48 GB | 48 h |

### Variant filtration

Variants are filtered using GATK `VariantFiltration` with the following criteria:

- **FS > 30.0** — Fisher strand bias
- **QD < 2.0** — Quality by depth
- **Cluster window:** 35 bp, cluster size 3

Variants failing these filters are flagged as `FSQD` in the `FILTER` column.

### Cleanup

Once all final outputs are verified (recalibrated BAM, filtered VCF, filtered TSV), intermediate files (read-group BAM, duplicate-marked BAM, split BAM, temp directories) are removed automatically to save scratch space.

---

## Output files

For each sample, the final outputs in `<SAMPLE>.processed/` are:

| File | Description |
|------|-------------|
| `<SAMPLE>.recal.bam` | Base quality score recalibrated BAM |
| `<SAMPLE>.raw.vcf` | Raw variants from HaplotypeCaller |
| `<SAMPLE>.filtered.vcf` | Variants after VariantFiltration |
| `<SAMPLE>.filtered.tsv` | Tab-separated variant table (CHROM, POS, TYPE, REF, ALT, QUAL, FILTER, allelic depths) |
| `<SAMPLE>.dup.metrics` | Picard MarkDuplicates metrics |
| `<SAMPLE>.recal.data.csv` | BQSR recalibration table |
---
