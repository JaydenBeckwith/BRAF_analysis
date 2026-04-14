#!/bin/bash
#
# RNA-Seq Realignment + Variant Calling Pipeline — PBS
#
# Three stages:
#   1) Build STAR index from chr* GRCh38 FASTA + GENCODE GTF (single job)
#   2) FASTQ → STAR two-pass alignment (per-sample jobs, depend on index)
#   3) Variant calling pipeline: AddRG → MarkDup → SplitNCigar → BQSR → HC → Filter
#      (per-sample jobs, depend on alignment)
#
# Usage:
#   qsub RNA_variant_pipeline.sh
#
# 
#
# *** CHECKPOINT SUPPORT ***
# Each step checks for its expected output file before running.
# If the output already exists (non-empty), the step is skipped.
# Timed-out or failed jobs can be resubmitted and will resume automatically.

set -euo pipefail

# ===================================================================
# CONFIGURATION
# ===================================================================
############
# Change your input paths and settings here before running. See comments for details.
############
# --- References ---
REF="/scratch/jo11/GRCh38.fa/GRCh38.primary_assembly.genome.fa"
GTF="/scratch/jo11/GRCh38.fa/gencode.v46.chr_patch_hapl_scaff.annotation.gtf"
STAR_INDEX_DIR="/scratch/jo11/STAR_index_GRCh38_gencode46"

# --- Input BAMs (original, pre-Picard) ---
BAM_FOLDER="/scratch/jo11/NeoTrio_RNA/CAGRF220510736_rnaseq"

# --- Singularity images ---
STAR_SIF="/scratch/jo11/singularity_images/star_2.7.11b--h43eeafb_0.sif"
PICARD_SIF="/scratch/jo11/singularity_images/picard.sif"
GATK_SIF="/scratch/jo11/singularity_images/gatk_4.5.0.0.sif"
SAMTOOLS_SIF="/scratch/jo11/singularity_images/samtools_v1.3.1_cv4.sif"

# --- Known sites for BQSR ---
MILLS="/scratch/jo11/NeoTrio_RNA/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KG="/scratch/jo11/NeoTrio_RNA/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP="/scratch/jo11/NeoTrio_RNA/dbsnp_138.hg38.vcf.gz"

# --- BRAF V600E locus (hg38) ---
BRAF_CHROM="chr7"
BRAF_POS=140753336

# --- Output directories ---
FASTQ_DIR="/scratch/jo11/NeoTrio_RNA_variants/fastqs"
ALIGN_DIR="/scratch/jo11/NeoTrio_RNA_variants/star_aligned"
VARIANT_DIR="/scratch/jo11/NeoTrio_RNA_variants"

# --- PBS settings ---
PBS_PROJECT="jo11"
PBS_QUEUE="normal"
PBS_STORAGE="scratch/jo11+gdata/jo11"
MAX_CONCURRENT=30

# --- Samples to skip (put exact BAM basenames without .bam, space-separated) ---
SKIP_SAMPLES=""

# ===================================================================
# SETUP
# ===================================================================

pbs_log_dir="${VARIANT_DIR}/pbs_logs"
mkdir -p "$pbs_log_dir" "$FASTQ_DIR" "$ALIGN_DIR"

bam_files=("$BAM_FOLDER"/*.bam)
bam_count=${#bam_files[@]}

echo ""
echo "================================================================="
echo -e "[$(date)]\tRNA-Seq Realignment + Variant Pipeline"
echo -e "  BAM folder:        $BAM_FOLDER"
echo -e "  Samples found:     $bam_count"
echo -e "  STAR index dir:    $STAR_INDEX_DIR"
echo -e "  Reference:         $REF"
echo -e "  GTF:               $GTF"
echo -e "  Skipping:          $SKIP_SAMPLES"
echo "================================================================="
echo ""

# ===================================================================
# STAGE 1: Build STAR index
# ===================================================================

echo "--- Stage 1: Submitting STAR index build ---"

index_job_id=$(qsub <<'INDEXEOF'
#!/bin/bash
#PBS -N star_index_build
#PBS -P jo11
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=24
#PBS -l mem=128GB
#PBS -l storage=scratch/jo11+gdata/jo11
#PBS -l wd
#PBS -o /scratch/jo11/NeoTrio_RNA_variants/pbs_logs/star_index.o
#PBS -e /scratch/jo11/NeoTrio_RNA_variants/pbs_logs/star_index.e

set -euo pipefail
module load singularity

STAR_SIF="/scratch/jo11/singularity_images/star_2.7.11b--h43eeafb_0.sif"
REF="/scratch/jo11/GRCh38.fa/GRCh38.primary_assembly.genome.fa"
GTF="/scratch/jo11/GRCh38.fa/gencode.v46.chr_patch_hapl_scaff.annotation.gtf"
STAR_INDEX_DIR="/scratch/jo11/STAR_index_GRCh38_gencode46"

mkdir -p "$STAR_INDEX_DIR"

# Check if index already exists
if [ -f "${STAR_INDEX_DIR}/SA" ] && [ -s "${STAR_INDEX_DIR}/SA" ]; then
    echo "[$(date)] STAR index already exists — skipping build"
    exit 0
fi

echo "[$(date)] Building STAR index..."

singularity exec "$STAR_SIF" STAR \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$REF" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 149 \
    --runThreadN 24

echo "[$(date)] STAR index build complete"
INDEXEOF
)

echo "  STAR index job: $index_job_id"

# ===================================================================
# STAGE 2 + 3: Per-sample — BAM/FASTQ --> STAR align --> Variant calling
# ===================================================================

echo ""
echo "--- Stage 2+3: Submitting per-sample realignment + variant calling ---"
echo ""

# Sort BAMs by timepoint for submission order
sort_by_timepoint() {
    for bam in "${bam_files[@]}"; do
        local name=$(basename "$bam")
        local order=5
        case "$name" in
            *PRE*)  order=1 ;;
            *ED1*)  order=2 ;;
            *ED2*)  order=3 ;;
            *CLND*) order=4 ;;
        esac
        echo "${order} ${bam}"
    done | sort -k1,1n -k2 | cut -d' ' -f2-
}

sorted_bams=()
while IFS= read -r bam; do
    sorted_bams+=("$bam")
done < <(sort_by_timepoint)

job_ids=()

for input_bam in "${sorted_bams[@]}"; do
    bn=$(basename "$input_bam" .bam)

    # Skip specific samples
    if [[ "$SKIP_SAMPLES" == *"$bn"* ]]; then
        echo "  Skipping: ${bn}"
        continue
    fi

    # Skip already completed samples
    if [ -f "${VARIANT_DIR}/${bn}.processed/${bn}.filtered.vcf" ] && \
       [ -f "${VARIANT_DIR}/${bn}.processed/${bn}.recal.bam" ]; then
        echo "  Already done: ${bn}"
        continue
    fi

    job_id=$(qsub -W depend=afterok:${index_job_id} <<EOF
#!/bin/bash
#PBS -N rna_${bn}
#PBS -P ${PBS_PROJECT}
#PBS -q ${PBS_QUEUE}
#PBS -l walltime=48:00:00
#PBS -l ncpus=12
#PBS -l mem=48GB
#PBS -l storage=${PBS_STORAGE}
#PBS -l wd
#PBS -o ${pbs_log_dir}/${bn}.o
#PBS -e ${pbs_log_dir}/${bn}.e

set -euo pipefail
module load singularity

# ---------------------------------------------------------------
# Checkpoint helper
# ---------------------------------------------------------------
checkpoint_run() {
    local step_label="\$1"
    local output_file="\$2"
    shift 2

    if [ -f "\$output_file" ] && [ -s "\$output_file" ]; then
        echo -e "[\$(date)]\t\${bn} — \${step_label} — SKIPPED (exists)"
        return 0
    fi

    echo -e "[\$(date)]\t\${bn} — \${step_label} — Running..."

    if "\$@"; then
        echo -e "[\$(date)]\t\${bn} — \${step_label} — Completed"
    else
        local rc=\$?
        echo -e "[\$(date)]\t\${bn} — \${step_label} — FAILED (exit \${rc})"
        rm -f "\$output_file"
        exit \$rc
    fi
}

# ---------------------------------------------------------------
# Paths
# ---------------------------------------------------------------
STAR_SIF="${STAR_SIF}"
PICARD_SIF="${PICARD_SIF}"
GATK_SIF="${GATK_SIF}"
SAMTOOLS_SIF="${SAMTOOLS_SIF}"
REF="${REF}"
MILLS="${MILLS}"
KG="${KG}"
DBSNP="${DBSNP}"
BRAF_CHROM="${BRAF_CHROM}"
BRAF_POS=${BRAF_POS}
STAR_INDEX_DIR="${STAR_INDEX_DIR}"

input_bam="${input_bam}"
bn="${bn}"
cpus=12

picard="singularity exec \$PICARD_SIF java -jar /usr/picard/picard.jar"
gatk="singularity exec \$GATK_SIF gatk"
samtools="singularity exec \$SAMTOOLS_SIF samtools"
star="singularity exec \$STAR_SIF STAR"

# Output directories
fqdir="${FASTQ_DIR}"
aligndir="${ALIGN_DIR}/\${bn}"
opdir="${VARIANT_DIR}/\${bn}.processed"
tmpdir="${VARIANT_DIR}/tmp/\${bn}"
mkdir -p "\$fqdir" "\$aligndir" "\$opdir" "\$tmpdir"

pre="\${opdir}/\${bn}"

# ===============================================================
# STAGE 2: BAM/FASTQ to STAR two-pass alignment
# ===============================================================

### or use fastq as input if already exists (checkpointed)
if [ -f "\${fqdir}/\${bn}_R1.fastq.gz" ] && [ -s "\${fqdir}/\${bn}_R1.fastq.gz" ]; then
    echo -e "[\$(date)]\t\${bn} — Step 1-2/11 - SKIPPED (FASTQs exist)"
else
    checkpoint_run "Step 1/12 - Sort by queryname" "\${fqdir}/\${bn}.namesorted.bam" \
        singularity exec \$SAMTOOLS_SIF samtools sort \
        -n \
        -o "\${fqdir}/\${bn}.namesorted.bam" \
        "\$input_bam"

    checkpoint_run "Step 2/12 - SamToFastq" "\${fqdir}/\${bn}_R1.fastq.gz" \
        \$picard SamToFastq \
        I="\${fqdir}/\${bn}.namesorted.bam" \
        F="\${fqdir}/\${bn}_R1.fastq.gz" \
        F2="\${fqdir}/\${bn}_R2.fastq.gz" \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR="\$tmpdir"

    rm -f "\${fqdir}/\${bn}.namesorted.bam"
fi

# ---------------------------------------------------------------
# Step 3/12 - STAR two-pass alignment
# ---------------------------------------------------------------
checkpoint_run "Step 3/12 - STAR two-pass alignment" "\${aligndir}/\${bn}.Aligned.sortedByCoord.out.bam" \
    \$star \
    --runThreadN \$cpus \
    --genomeDir "\$STAR_INDEX_DIR" \
    --readFilesIn "\${fqdir}/\${bn}_R1.fastq.gz" "\${fqdir}/\${bn}_R2.fastq.gz" \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 40000000000 \
    --outFileNamePrefix "\${aligndir}/\${bn}." \
    --outTmpDir "\${tmpdir}/\${bn}_STARtmp" \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 2

# Index the STAR BAM
if [ ! -f "\${aligndir}/\${bn}.Aligned.sortedByCoord.out.bam.bai" ]; then
    singularity exec \$SAMTOOLS_SIF samtools index \
        "\${aligndir}/\${bn}.Aligned.sortedByCoord.out.bam"
fi

# Use the new STAR-aligned BAM for downstream steps
star_bam="\${aligndir}/\${bn}.Aligned.sortedByCoord.out.bam"

# ===============================================================
# STAGE 3: Variant calling pipeline
# ===============================================================
# ---------------------------------------------------------------
# Step 4 — AddOrReplaceReadGroups (Picard)
# ---------------------------------------------------------------
checkpoint_run "Step 4/12 - AddOrReplaceReadGroups" "\${pre}.rg.bam" \
    \$picard AddOrReplaceReadGroups \
    I="\$star_bam" \
    O="\${pre}.rg.bam" \
    RGID="\${bn}" \
    RGLB="lib1" \
    RGPL="ILLUMINA" \
    RGPU="unit1" \
    RGSM="\${bn}" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="\$tmpdir"

# ---------------------------------------------------------------
# Step 5 — Mark Duplicates (Picard)
# ---------------------------------------------------------------
checkpoint_run "Step 5/12 - MarkDuplicates" "\${pre}.dupMarked.bam" \
    \$picard MarkDuplicates \
    I="\${pre}.rg.bam" \
    O="\${pre}.dupMarked.bam" \
    M="\${pre}.dup.metrics" \
    TMP_DIR="\$tmpdir" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

# ---------------------------------------------------------------
# Step 6 — SplitNCigarReads (GATK)
# ---------------------------------------------------------------
checkpoint_run "Step 6/12 - SplitNCigarReads" "\${pre}.split.bam" \
    \$gatk SplitNCigarReads \
    -R "\$REF" \
    -I "\${pre}.dupMarked.bam" \
    -O "\${pre}.split.bam" \
    --tmp-dir "\$tmpdir"

# ---------------------------------------------------------------
# Step 7 — BaseRecalibrator
# ---------------------------------------------------------------
checkpoint_run "Step 7/12 - BaseRecalibrator" "\${pre}.recal.data.csv" \
    \$gatk BaseRecalibrator \
    -R "\$REF" \
    -I "\${pre}.split.bam" \
    -O "\${pre}.recal.data.csv" \
    --known-sites "\$KG" \
    --known-sites "\$MILLS" \
    --known-sites "\$DBSNP"

# ---------------------------------------------------------------
# Step 8 — ApplyBQSR
# ---------------------------------------------------------------
checkpoint_run "Step 8/12 - ApplyBQSR" "\${pre}.recal.bam" \
    \$gatk ApplyBQSR \
    -R "\$REF" \
    -I "\${pre}.split.bam" \
    -O "\${pre}.recal.bam" \
    --bqsr-recal-file "\${pre}.recal.data.csv"

# ---------------------------------------------------------------
# Step 9 — HaplotypeCaller
# ---------------------------------------------------------------
checkpoint_run "Step 9/12 - HaplotypeCaller" "\${pre}.raw.vcf" \
    \$gatk HaplotypeCaller \
    -R "\$REF" \
    -I "\${pre}.recal.bam" \
    -O "\${pre}.raw.vcf" \
    -D "\$DBSNP" \
    -stand-call-conf 20.0 \
    --dont-use-soft-clipped-bases \
    --native-pair-hmm-threads \$cpus

# ---------------------------------------------------------------
# Step 10 — VariantFiltration
# ---------------------------------------------------------------
checkpoint_run "Step 10/12 - VariantFiltration" "\${pre}.filtered.vcf" \
    \$gatk VariantFiltration \
    -R "\$REF" \
    -V "\${pre}.raw.vcf" \
    -O "\${pre}.filtered.vcf" \
    -window 35 \
    -cluster 3 \
    -filter "FS > 30.0 || QD < 2.0" \
    --filter-name FSQD

# ---------------------------------------------------------------
# Step 11 — VariantsToTable
# ---------------------------------------------------------------
checkpoint_run "Step 11 - VariantsToTable" "\${pre}.filtered.tsv" \
    \$gatk VariantsToTable \
    -V "\${pre}.filtered.vcf" \
    -F CHROM -F POS -F TYPE -F REF -F ALT -F QUAL -F FILTER -GF AD \
    -O "\${pre}.filtered.tsv"


# ---------------------------------------------------------------
# Cleanup intermediates to save space
# ---------------------------------------------------------------
all_done=true
for f in "\${pre}.recal.bam" "\${pre}.filtered.vcf" "\${pre}.filtered.tsv" "\${pre}.braf_allelic_counts.tsv"; do
    if [ ! -f "\$f" ] || [ ! -s "\$f" ]; then
        all_done=false
        break
    fi
done

if [ "\$all_done" = true ]; then
    rm -f "\${pre}.rg.bam" "\${pre}.rg.bai"
    rm -f "\${pre}.dupMarked.bam" "\${pre}.dupMarked.bai"
    rm -f "\${pre}.split.bam" "\${pre}.split.bai"
    rm -f "\${fqdir}/\${bn}.namesorted.bam"
    rm -rf "\$tmpdir"
    echo -e "[\$(date)]\t\${bn} — DONE (all outputs verified, intermediates cleaned)"
else
    echo -e "[\$(date)]\t\${bn} — WARNING: some outputs missing, intermediates preserved"
fi
EOF
)

    job_ids+=("$job_id")
    echo "  Submitted: ${bn}  →  ${job_id}"

    # Throttle
    if [ "${#job_ids[@]}" -ge "$MAX_CONCURRENT" ]; then
        oldest="${job_ids[0]}"
        echo "  (Reached $MAX_CONCURRENT concurrent — waiting on ${oldest})"
        while qstat "$oldest" &>/dev/null 2>&1; do
            sleep 60
        done
        job_ids=("${job_ids[@]:1}")
    fi
done

echo ""
echo "================================================================="
echo "  All $bam_count sample jobs submitted (depend on index job: $index_job_id)"
echo "================================================================="

echo ""
echo "================================================================="
echo -e "[\$(date)]\tPIPELINE SUMMARY"
echo -e "  Total samples:     \$total"
echo -e "  Completed:         \$completed"
echo -e "  Failed/incomplete: \$((total - completed))"
echo ""

if [ -n "\$failed_samples" ]; then
    echo "  Failed samples:"
    echo -e "\$failed_samples" | sed 's/^/    /'
fi

echo ""
echo "================================================================="
echo ""
COLLECTOREOF
