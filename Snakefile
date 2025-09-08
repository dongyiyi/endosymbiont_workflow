# Snakemake workflow: QC → Cutadapt → QIIME2 (DADA2) → Taxonomy (SILVA + GG2) → Concordance
# Run:
#   snakemake -j 8 --use-conda --configfile workflow/config.yaml
# ================================================================
# Endosymbiont 16S Amplicon Workflow - Snakefile (full, updated)
# ================================================================

import os, re, sys
from pathlib import Path

# ---------------- Config / Paths ----------------
cfg = config

RAW_DIR = cfg.get("raw_dir", "raw")
META    = Path(cfg.get("metadata_tsv", "metadata/sample-metadata.tsv"))

R1SFX = cfg.get("r1_suffix", "_R1.fastq.gz")
R2SFX = cfg.get("r2_suffix", "_R2.fastq.gz")

CUTADAPT_EXTRA = cfg.get("cutadapt_extra", "--minimum-length 50 -q 20,20")

# DADA2 params
dada2_cfg = cfg.get("dada2", {})
TRIM_LEFT_F  = int(dada2_cfg.get("trim_left_f", 0))
TRIM_LEFT_R  = int(dada2_cfg.get("trim_left_r", 0))
TRUNC_LEN_F  = int(dada2_cfg.get("trunc_len_f", 0))
TRUNC_LEN_R  = int(dada2_cfg.get("trunc_len_r", 0))

# QIIME env selection
USE_EXISTING_QIIME = bool(cfg.get("use_existing_qiime_env", True))
QIIME_ENV_NAME     = cfg.get("qiime_env_name", "qiime2-amplicon-2024.10")
Q2 = ('conda run -n ' + QIIME_ENV_NAME + ' qiime') if USE_EXISTING_QIIME else 'qiime'

# Primary (SILVA) classifier
PRIMARY_CLASSIFIER = Path(cfg.get("classifier_qza", "refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-classifier.qza"))

# Reference build params (used by build_classifier_silva)
ref_cfg = cfg.get("reference", {})
region_label  = ref_cfg.get("region_label", "V3V4")
silva_version = str(ref_cfg.get("silva_version", "138.2"))
silva_target  = ref_cfg.get("silva_target", "SSURef_NR99")
primers       = ref_cfg.get("primers", {})
prim_fwd      = primers.get("forward", "CCTACGGGNGGCWGCAG")
prim_rev      = primers.get("reverse", "GACTACHVGGGTATCTAATCC")
symbiont_list = ref_cfg.get("symbiont_list", "refs/symbiont_taxa.txt")

# GG2 (optional)
GG2 = cfg.get("gg2", {})
GG2_DIR  = Path(GG2.get("outdir", "refs/GG2-V3V4"))
GG2_SEQS = GG2.get("seqs_qza")           # refs/GG2-raw/2024.09.backbone.full-length.fna.qza (used to build classifier)
GG2_TAX  = GG2.get("tax_qza")            # refs/GG2-raw/2024.09.backbone.tax.qza
GG2_SEPP_REF = GG2.get("sepp_ref_qza")   # refs/GG2-raw/2022.10.backbone.sepp-reference.qza  (optional)
use_gg2 = bool(GG2_SEQS) and bool(GG2_TAX)

# For OTU-97 with GG2 (closed/open reference) we need the  V3V4-uniq-seqs.qza built by your gg2 classifier rule
if use_gg2:
    GG2_REF_SEQS = str(GG2_DIR / "V3V4-uniq-seqs.qza")

# Summary options
summary_cfg = cfg.get("summary", {})
GROUP_COL   = summary_cfg.get("group_column", "group")

# Concordance options
TOP_N = int(cfg.get("concordance", {}).get("top_n", 50))

# --------------- Helper: sample list ---------------
def _strip_suffix(name, sfx):
    return name[:-len(sfx)] if name.endswith(sfx) else None

r1_bases = { _strip_suffix(p.name, R1SFX) for p in Path(RAW_DIR).glob(f"*{R1SFX}") }
r2_bases = { _strip_suffix(p.name, R2SFX) for p in Path(RAW_DIR).glob(f"*{R2SFX}") }
SAMPLES = sorted({s for s in r1_bases.intersection(r2_bases) if s})

if len(SAMPLES) == 0:
    # 
    raise AssertionError(f"No samples found in {RAW_DIR} matching *{R1SFX} and *{R2SFX}")

# --------------- Wildcard constraints ---------------
# (ASV/OTU )
wildcard_constraints:
    rank="phylum|family|genus|species"

# ====================================================
# QC (FastQC/MultiQC)  ——
# ====================================================
rule fastqc:
    input:
        r1=lambda wc: f"{RAW_DIR}/{wc.sample}{R1SFX}",
        r2=lambda wc: f"{RAW_DIR}/{wc.sample}{R2SFX}"
    output:
        html1=f"qiime2/qc/fastqc/{{sample}}_R1_fastqc.html",
        html2=f"qiime2/qc/fastqc/{{sample}}_R2_fastqc.html",
        zip1 =f"qiime2/qc/fastqc/{{sample}}_R1_fastqc.zip",
        zip2 =f"qiime2/qc/fastqc/{{sample}}_R2_fastqc.zip"
    conda: "workflow/envs/qc.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p qiime2/qc/fastqc
        fastqc -t {threads} -o qiime2/qc/fastqc {input.r1} {input.r2}
        """

rule multiqc:
    input:
        expand("qiime2/qc/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("qiime2/qc/fastqc/{sample}_R2_fastqc.zip", sample=SAMPLES)
    output:
        html="qiime2/qc/multiqc/multiqc_report.html"
    conda: "workflow/envs/qc.yaml"
    shell:
        r"""
        mkdir -p qiime2/qc/multiqc
        multiqc -o qiime2/qc/multiqc qiime2/qc/fastqc
        """

# ====================================================
# Trim (cutadapt) → Manifest → QIIME Import → DADA2
# ====================================================
rule cutadapt_paired:
    input:
        r1=lambda wc: f"{RAW_DIR}/{wc.sample}{R1SFX}",
        r2=lambda wc: f"{RAW_DIR}/{wc.sample}{R2SFX}"
    output:
        r1=f"qiime2/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"qiime2/trimmed/{{sample}}_R2.trimmed.fastq.gz"
    conda: "workflow/envs/cutadapt.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p qiime2/trimmed
        cutadapt -j {threads} \
          -g {prim_fwd} -G {prim_rev} \
          -o {output.r1} -p {output.r2} \
          {CUTADAPT_EXTRA} \
          {input.r1} {input.r2} > qiime2/trimmed/{wildcards.sample}.cutadapt.log
        """

rule make_manifest:
    input:
        r1=expand("qiime2/trimmed/{sample}_R1.trimmed.fastq.gz", sample=SAMPLES),
        r2=expand("qiime2/trimmed/{sample}_R2.trimmed.fastq.gz", sample=SAMPLES)
    output:
        "qiime2/manifest_pe.tsv"
    run:
        Path("qiime2").mkdir(exist_ok=True, parents=True)
        with open("qiime2/manifest_pe.tsv", "w") as f:
            print("sample-id\tforward-absolute-filepath\treverse-absolute-filepath", file=f)
            for s in SAMPLES:
                fwd = Path(f"qiime2/trimmed/{s}_R1.trimmed.fastq.gz").resolve()
                rev = Path(f"qiime2/trimmed/{s}_R2.trimmed.fastq.gz").resolve()
                print(f"{s}\t{fwd}\t{rev}", file=f)

rule qiime_import:
    input: "qiime2/manifest_pe.tsv"
    output: "qiime2/demux.qza"
    shell:
        r"""
        {Q2} tools import \
          --type 'SampleData[PairedEndSequencesWithQuality]' \
          --input-format PairedEndFastqManifestPhred33V2 \
          --input-path {input} \
          --output-path {output}
        """

rule demux_summarize:
    input: "qiime2/demux.qza"
    output: "qiime2/demux.qzv"
    shell:
        r"""
        {Q2} demux summarize --i-data {input} --o-visualization {output}
        """

rule dada2_denoise_paired:
    input: "qiime2/demux.qza"
    output:
        table="qiime2/table.qza",
        reps ="qiime2/rep-seqs.qza",
        stats="qiime2/denoise-stats.qza"
    threads: 8
    shell:
        r"""
        {Q2} dada2 denoise-paired \
          --i-demultiplexed-seqs {input} \
          --p-trim-left-f {TRIM_LEFT_F} --p-trim-left-r {TRIM_LEFT_R} \
          --p-trunc-len-f {TRUNC_LEN_F} --p-trunc-len-r {TRUNC_LEN_R} \
          --o-table {output.table} \
          --o-representative-sequences {output.reps} \
          --o-denoising-stats {output.stats} \
          --p-n-threads {threads}
        """

# ====================================================
# Taxonomy (Primary SILVA) + Barplot
# ====================================================
rule taxonomy_primary:
    input:
        reps="qiime2/rep-seqs.qza",
        classifier=str(PRIMARY_CLASSIFIER),
        meta=str(META)
    output:
        tax="qiime2/taxonomy/taxonomy_primary.qza",
        bar="qiime2/taxonomy/taxa-barplot_primary.qzv"
    shell:
        r"""
        mkdir -p qiime2/taxonomy
        {Q2} feature-classifier classify-sklearn \
          --i-classifier {input.classifier} \
          --i-reads {input.reps} \
          --o-classification {output.tax}
        {Q2} taxa barplot \
          --i-table qiime2/table.qza \
          --i-taxonomy {output.tax} \
          --m-metadata-file "{input.meta}" \
          --o-visualization {output.bar}
        """

# ----------------------------------------------------
# TSV Exports (always available)
# ----------------------------------------------------
rule export_taxonomy_primary_tsv:
    input: "qiime2/taxonomy/taxonomy_primary.qza"
    output: tsv="qiime2/taxonomy/taxonomy_primary.tsv"
    shell:
        r"""
        rm -rf qiime2/taxonomy/_export_primary || true
        {Q2} tools export --input-path {input} --output-path qiime2/taxonomy/_export_primary
        mv qiime2/taxonomy/_export_primary/taxonomy.tsv {output.tsv}
        rm -rf qiime2/taxonomy/_export_primary
        """

rule export_table_tsv:
    input: "qiime2/table.qza"
    output:
        biom="qiime2/_tmp/feature-table.biom",
        tsv ="qiime2/feature-table.tsv"
    shell:
        r"""
        rm -rf qiime2/_tmp || true; mkdir -p qiime2/_tmp
        {Q2} tools export --input-path {input} --output-path qiime2/_tmp
        conda run -n {QIIME_ENV_NAME} biom convert -i {output.biom} -o {output.tsv} --to-tsv
        """

# ====================================================
# ASV 
# ====================================================
rule asv_level_table:
    input:
        table_qza="qiime2/table.qza",
        tax_qza="qiime2/taxonomy/taxonomy_primary.qza"
    output:
        tsv = "qiime2/summary/asv_{rank}.tsv"
    params:
        rank=lambda wc: wc.rank
    shell:
        r"""
        set -euo pipefail
        TMP=qiime2/summary/_tmp_asv_{wildcards.rank}
        mkdir -p "$TMP/tab" "$TMP/tax"
        {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
        conda run -n {QIIME_ENV_NAME} biom convert -i "$TMP/tab/feature-table.biom" -o "$TMP/feature-table.tsv" --to-tsv
        {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

        conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
          --table-tsv "$TMP/feature-table.tsv" \
          --tax-tsv   "$TMP/tax/taxonomy.tsv" \
          --rank {params.rank} \
          --output {output.tsv}
        """

rule asv_post_rank:
    input:
        wide = "qiime2/summary/asv_{rank}.tsv",
        meta = str(META)
    output:
        rel  = "qiime2/summary/asv_{rank}.relabund.tsv",
        long = "qiime2/summary/asv_{rank}.long.tsv"
    params:
        group = GROUP_COL
    shell:
        r"""
        conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
          --wide-tsv {input.wide} \
          --metadata {input.meta} \
          --group-col {params.group} \
          --out-rel  {output.rel} \
          --out-long {output.long}
        """

# ====================================================
# OTU97（de novo, vsearch）→ Primary SILVA 
# ====================================================
rule otu97_cluster_denovo:
    input:
        table="qiime2/table.qza",
        reps ="qiime2/rep-seqs.qza"
    output:
        table="qiime2/otu97/otu97-table.qza",
        reps ="qiime2/otu97/otu97-rep-seqs.qza"
    shell:
        r"""
        mkdir -p qiime2/otu97
        {Q2} vsearch cluster-features-de-novo \
          --i-table {input.table} \
          --i-sequences {input.reps} \
          --p-perc-identity 0.97 \
          --o-clustered-table {output.table} \
          --o-clustered-sequences {output.reps}
        """

rule otu97_taxonomy_primary:
    input:
        reps="qiime2/otu97/otu97-rep-seqs.qza",
        classifier=str(PRIMARY_CLASSIFIER)
    output:
        tax="qiime2/otu97/taxonomy_otu97_primary.qza"
    shell:
        r"""
        {Q2} feature-classifier classify-sklearn \
          --i-classifier {input.classifier} \
          --i-reads {input.reps} \
          --o-classification {output.tax}
        """

rule otu97_primary_level_table:
    input:
        table_qza="qiime2/otu97/otu97-table.qza",
        tax_qza  ="qiime2/otu97/taxonomy_otu97_primary.qza"
    output:
        tsv="qiime2/summary/otu97_{rank}.tsv"
    params:
        rank=lambda wc: wc.rank
    shell:
        r"""
        set -euo pipefail
        TMP=qiime2/summary/_tmp_otu97_{wildcards.rank}
        mkdir -p "$TMP/tab" "$TMP/tax"
        {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
        conda run -n {QIIME_ENV_NAME} biom convert -i "$TMP/tab/feature-table.biom" -o "$TMP/otu-table.tsv" --to-tsv
        {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

        conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
          --table-tsv "$TMP/otu-table.tsv" \
          --tax-tsv   "$TMP/tax/taxonomy.tsv" \
          --rank {params.rank} \
          --output {output.tsv}
        """

rule otu97_primary_post_rank:
    input:
        wide = "qiime2/summary/otu97_{rank}.tsv",
        meta = str(META)
    output:
        rel  = "qiime2/summary/otu97_{rank}.relabund.tsv",
        long = "qiime2/summary/otu97_{rank}.long.tsv"
    params:
        group = GROUP_COL
    shell:
        r"""
        conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
          --wide-tsv {input.wide} \
          --metadata {input.meta} \
          --group-col {params.group} \
          --out-rel  {output.rel} \
          --out-long {output.long}
        """

# ====================================================
# GG2 Cross-check（ASV/OTU97 with GG2 taxonomy）
# ====================================================
if use_gg2:

    rule taxonomy_gg2:
        input:
            reps="qiime2/rep-seqs.qza",
            classifier=str(GG2_DIR / "V3V4-uniq-classifier.qza")
        output:
            tax="qiime2/taxonomy/taxonomy_gg2.qza"
        shell:
            r"""
            mkdir -p qiime2/taxonomy
            {Q2} feature-classifier classify-sklearn \
              --i-classifier {input.classifier} \
              --i-reads {input.reps} \
              --o-classification {output.tax}
            """

    rule export_taxonomy_gg2_tsv:
        input: "qiime2/taxonomy/taxonomy_gg2.qza"
        output: tsv="qiime2/taxonomy/taxonomy_gg2.tsv"
        shell:
            r"""
            rm -rf qiime2/taxonomy/_export_gg2 || true
            {Q2} tools export --input-path {input} --output-path qiime2/taxonomy/_export_gg2
            mv qiime2/taxonomy/_export_gg2/taxonomy.tsv {output.tsv}
            rm -rf qiime2/taxonomy/_export_gg2
            """

    # ---- ASV with GG2 taxonomy ----
    rule asv_level_table_gg2:
        input:
            table_qza="qiime2/table.qza",
            tax_qza  ="qiime2/taxonomy/taxonomy_gg2.qza"
        output:
            tsv = "qiime2/summary/asv_gg2_{rank}.tsv"
        params:
            rank=lambda wc: wc.rank
        shell:
            r"""
            set -euo pipefail
            TMP=qiime2/summary/_tmp_asv_gg2_{wildcards.rank}
            mkdir -p "$TMP/tab" "$TMP/tax"
            {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
            conda run -n {QIIME_ENV_NAME} biom convert -i "$TMP/tab/feature-table.biom" -o "$TMP/feature-table.tsv" --to-tsv
            {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

            conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
              --table-tsv "$TMP/feature-table.tsv" \
              --tax-tsv   "$TMP/tax/taxonomy.tsv" \
              --rank {params.rank} \
              --output {output.tsv}
            """

    rule asv_gg2_post_rank:
        input:
            wide = "qiime2/summary/asv_gg2_{rank}.tsv",
            meta = str(META)
        output:
            rel  = "qiime2/summary/asv_gg2_{rank}.relabund.tsv",
            long = "qiime2/summary/asv_gg2_{rank}.long.tsv"
        params:
            group = GROUP_COL
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
              --wide-tsv {input.wide} \
              --metadata {input.meta} \
              --group-col {params.group} \
              --out-rel  {output.rel} \
              --out-long {output.long}
            """

    # ---- OTU97 using GG2 reference (closed/open) ----
    # closed-reference
    rule otu97gg2_cluster_closed:
        input:
            table="qiime2/table.qza",
            reps ="qiime2/rep-seqs.qza",
            ref  = GG2_REF_SEQS
        output:
            table="qiime2/otu97gg2/closed-table.qza",
            reps ="qiime2/otu97gg2/closed-rep-seqs.qza"
        shell:
            r"""
            mkdir -p qiime2/otu97gg2
            {Q2} vsearch cluster-features-closed-reference \
              --i-table {input.table} \
              --i-sequences {input.reps} \
              --i-reference-sequences {input.ref} \
              --p-perc-identity 0.97 \
              --o-clustered-table {output.table} \
              --o-clustered-sequences {output.reps}
            """

    # open-reference
    rule otu97gg2_cluster_open:
        input:
            table="qiime2/table.qza",
            reps ="qiime2/rep-seqs.qza",
            ref  = GG2_REF_SEQS
        output:
            table="qiime2/otu97gg2/open-table.qza",
            reps ="qiime2/otu97gg2/open-rep-seqs.qza"
        shell:
            r"""
            mkdir -p qiime2/otu97gg2
            {Q2} vsearch cluster-features-open-reference \
              --i-table {input.table} \
              --i-sequences {input.reps} \
              --i-reference-sequences {input.ref} \
              --p-perc-identity 0.97 \
              --o-clustered-table {output.table} \
              --o-clustered-sequences {output.reps}
            """

    # GG2 taxonomy on closed/open OTU reps
    rule otu97gg2_taxonomy_closed:
        input:
            reps="qiime2/otu97gg2/closed-rep-seqs.qza",
            classifier=str(GG2_DIR / "V3V4-uniq-classifier.qza")
        output:
            tax="qiime2/otu97gg2/taxonomy_closed_gg2.qza"
        shell:
            r"""
            {Q2} feature-classifier classify-sklearn \
              --i-classifier {input.classifier} \
              --i-reads {input.reps} \
              --o-classification {output.tax}
            """

    rule otu97gg2_taxonomy_open:
        input:
            reps="qiime2/otu97gg2/open-rep-seqs.qza",
            classifier=str(GG2_DIR / "V3V4-uniq-classifier.qza")
        output:
            tax="qiime2/otu97gg2/taxonomy_open_gg2.qza"
        shell:
            r"""
            {Q2} feature-classifier classify-sklearn \
              --i-classifier {input.classifier} \
              --i-reads {input.reps} \
              --o-classification {output.tax}
            """

    # OTU97 (GG2)
    rule otu97_gg2_level_table_closed:
        input:
            table_qza="qiime2/otu97gg2/closed-table.qza",
            tax_qza  ="qiime2/otu97gg2/taxonomy_closed_gg2.qza"
        output:
            tsv="qiime2/summary/otu97_gg2_closed_{rank}.tsv"
        params:
            rank=lambda wc: wc.rank
        shell:
            r"""
            set -euo pipefail
            TMP=qiime2/summary/_tmp_otu97_gg2_closed_{wildcards.rank}
            mkdir -p "$TMP/tab" "$TMP/tax"
            {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
            conda run -n {QIIME_ENV_NAME} biom convert -i "$TMP/tab/feature-table.biom" -o "$TMP/otu-table.tsv" --to-tsv
            {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

            conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
              --table-tsv "$TMP/otu-table.tsv" \
              --tax-tsv   "$TMP/tax/taxonomy.tsv" \
              --rank {params.rank} \
              --output {output.tsv}
            """

    rule otu97_gg2_level_table_open:
        input:
            table_qza="qiime2/otu97gg2/open-table.qza",
            tax_qza  ="qiime2/otu97gg2/taxonomy_open_gg2.qza"
        output:
            tsv="qiime2/summary/otu97_gg2_open_{rank}.tsv"
        params:
            rank=lambda wc: wc.rank
        shell:
            r"""
            set -euo pipefail
            TMP=qiime2/summary/_tmp_otu97_gg2_open_{wildcards.rank}
            mkdir -p "$TMP/tab" "$TMP/tax"
            {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
            conda run -n {QIIME_ENV_NAME} biom convert -i "$TMP/tab/feature-table.biom" -o "$TMP/otu-table.tsv" --to-tsv
            {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

            conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
              --table-tsv "$TMP/otu-table.tsv" \
              --tax-tsv   "$TMP/tax/taxonomy.tsv" \
              --rank {params.rank} \
              --output {output.tsv}
            """

    rule otu97_gg2_post_rank_closed:
        input:
            wide = "qiime2/summary/otu97_gg2_closed_{rank}.tsv",
            meta = str(META)
        output:
            rel  = "qiime2/summary/otu97_gg2_closed_{rank}.relabund.tsv",
            long = "qiime2/summary/otu97_gg2_closed_{rank}.long.tsv"
        params:
            group = GROUP_COL
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
              --wide-tsv {input.wide} \
              --metadata {input.meta} \
              --group-col {params.group} \
              --out-rel  {output.rel} \
              --out-long {output.long}
            """

    rule otu97_gg2_post_rank_open:
        input:
            wide = "qiime2/summary/otu97_gg2_open_{rank}.tsv",
            meta = str(META)
        output:
            rel  = "qiime2/summary/otu97_gg2_open_{rank}.relabund.tsv",
            long = "qiime2/summary/otu97_gg2_open_{rank}.long.tsv"
        params:
            group = GROUP_COL
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
              --wide-tsv {input.wide} \
              --metadata {input.meta} \
              --group-col {params.group} \
              --out-rel  {output.rel} \
              --out-long {output.long}
            """

    # Concordance top-N（taxonomy TSV）
    rule concordance_top:
        input:
            table_tsv="qiime2/feature-table.tsv",
            silva_tsv="qiime2/taxonomy/taxonomy_primary.tsv",
            gg2_tsv  ="qiime2/taxonomy/taxonomy_gg2.tsv"
        output:
            csv=f"qiime2/concordance/top{TOP_N}_concordance.csv"
        conda: "workflow/envs/qc.yaml"
        script:
            "workflow/scripts/concordance.py"

    # SEPP insertion tree with GG2 non-V4 backbone
    if GG2_SEPP_REF:
        rule gg2_nonv4_sepp_insertion:
            input:
                reps="qiime2/rep-seqs.qza",
                ref =GG2_SEPP_REF
            output:
                tree="qiime2/trees/gg2_nonv4_sepp/insertion-tree.qza",
                place="qiime2/trees/gg2_nonv4_sepp/placements.qza"
            shell:
                r"""
                mkdir -p qiime2/trees/gg2_nonv4_sepp
                {Q2} fragment-insertion sepp \
                  --i-representative-sequences {input.reps} \
                  --i-reference-database {input.ref} \
                  --o-tree {output.tree} \
                  --o-placements {output.place}
                """

        rule export_newick_gg2_sepp:
            input: "qiime2/trees/gg2_nonv4_sepp/insertion-tree.qza"
            output: "qiime2/trees/gg2_nonv4_sepp/tree.nwk"
            shell:
                r"""
                rm -rf qiime2/trees/gg2_nonv4_sepp/_exp || true
                {Q2} tools export --input-path {input} --output-path qiime2/trees/gg2_nonv4_sepp/_exp
                mv qiime2/trees/gg2_nonv4_sepp/_exp/tree.nwk {output}
                rm -rf qiime2/trees/gg2_nonv4_sepp/_exp
                """

# ====================================================
# （SILVA / GG2）
# ====================================================
# —— SILVA ——  refs/get_silva_classifier.sh 
rule build_classifier_silva:
    output:
        classifier=f"refs/{region_label}-silva-{silva_version}/silva-{silva_version}-{region_label}-uniq-classifier.qza",
        uniq_seqs =f"refs/{region_label}-silva-{silva_version}/silva-{silva_version}-{region_label}-uniq-seqs.qza",
        uniq_tax  =f"refs/{region_label}-silva-{silva_version}/silva-{silva_version}-{region_label}-uniq-tax.qza"
    params:
        region=region_label, ver=silva_version, target=silva_target,
        fwd=prim_fwd, rev=prim_rev, symb=symbiont_list
    shell:
        r"""
        TARGET="{params.target}" \
        bash refs/get_silva_classifier.sh {params.region} '{params.fwd}' '{params.rev}' {params.ver} {params.symb}
        """

# —— GG2 ——  refs/GG2-raw/*.qza 
rule build_classifier_gg2:
    input:
        seqs = str(GG2_SEQS),
        tax  = str(GG2_TAX)
    output:
        classifier=str(GG2_DIR / "V3V4-uniq-classifier.qza")
    params:
        region=region_label, outdir=str(GG2_DIR), fwd=prim_fwd, rev=prim_rev
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        QIIME="{Q2}"

        $QIIME feature-classifier extract-reads \
          --i-sequences "{input.seqs}" \
          --p-f-primer '{params.fwd}' \
          --p-r-primer '{params.rev}' \
          --o-reads "{params.outdir}/reads.qza"

        # dereplicate
        DEREP_FLAG="--o-dereplicated-taxa"
        $QIIME rescript dereplicate --help | grep -q -- "--o-dereplicated-taxonomy" && DEREP_FLAG="--o-dereplicated-taxonomy"

        $QIIME rescript dereplicate \
          --i-sequences "{params.outdir}/reads.qza" \
          --i-taxa      "{input.tax}" \
          --p-mode uniq \
          --o-dereplicated-sequences "{params.outdir}/V3V4-uniq-seqs.qza" \
          $DEREP_FLAG                       "{params.outdir}/V3V4-uniq-tax.qza"

        $QIIME feature-classifier fit-classifier-naive-bayes \
          --i-reference-reads    "{params.outdir}/V3V4-uniq-seqs.qza" \
          --i-reference-taxonomy "{params.outdir}/V3V4-uniq-tax.qza" \
          --o-classifier         "{output.classifier}"
        """

# ====================================================
# Convenience: rule all (optional)
# ====================================================
rule all:
    input:
        # 
        "qiime2/table.qza",
        "qiime2/rep-seqs.qza",
        "qiime2/demux.qzv",
        "qiime2/taxonomy/taxonomy_primary.qza",
        "qiime2/taxonomy/taxa-barplot_primary.qzv",
        # ASV 
        expand("qiime2/summary/asv_{rank}.tsv", rank=["phylum","family","genus","species"]),
        # ASV 
        expand("qiime2/summary/asv_{rank}.relabund.tsv", rank=["phylum","family","genus","species"]),
        expand("qiime2/summary/asv_{rank}.long.tsv",     rank=["phylum","family","genus","species"]),
        # OTU97 (de novo + primary)
        "qiime2/otu97/otu97-table.qza",
        "qiime2/otu97/otu97-rep-seqs.qza",
        "qiime2/otu97/taxonomy_otu97_primary.qza",
        expand("qiime2/summary/otu97_{rank}.tsv",         rank=["phylum","family","genus","species"]),
        expand("qiime2/summary/otu97_{rank}.relabund.tsv",rank=["phylum","family","genus","species"]),
        expand("qiime2/summary/otu97_{rank}.long.tsv",    rank=["phylum","family","genus","species"])
        
