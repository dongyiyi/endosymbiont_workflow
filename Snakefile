# Snakemake workflow: QC → Cutadapt → QIIME2 (DADA2) → Taxonomy (SILVA + GG2) → Concordance
# Run:
#   snakemake -j 8 --use-conda --configfile workflow/config.yaml

import yaml
from pathlib import Path

configfile: "workflow/config.yaml"
cfg = config

# 只允许 rank 取四个合法值，防止 asv_level_table 误匹配到 *.relabund.tsv / *.long.tsv
wildcard_constraints:
    rank="phylum|family|genus|species"

RAW = Path(cfg["raw_dir"]).resolve()
META = Path(cfg["metadata_tsv"]).resolve()
OUT = Path(cfg.get("outdir", "qiime2")).resolve()

R1SFX = cfg.get("r1_suffix", "_R1.fastq.gz")
R2SFX = cfg.get("r2_suffix", "_R2.fastq.gz")

SAMPLES = sorted([p.name.replace(R1SFX,"") for p in RAW.glob(f"*{R1SFX}")])
assert SAMPLES, f"No samples found in {RAW} matching *{R1SFX}"

GG2_SEPP_REF = cfg.get("gg2", {}).get("sepp_ref_qza", "")

USE_EXISTING_QIIME = bool(cfg.get("use_existing_qiime_env", False))
QIIME_ENV_NAME = cfg.get("qiime_env_name","qiime2-amplicon-2024.10")

# Convenience: run qiime consistently
Q2 = f"conda run -n {QIIME_ENV_NAME} qiime" if USE_EXISTING_QIIME else "qiime"

def qiime(cmd):
    if USE_EXISTING_QIIME:
        return f"conda run -n {QIIME_ENV_NAME} qiime {cmd}"
    else:
        return f"qiime {cmd}"

# --- Reference / classifier settings ---
ref = cfg.get("reference", {})
region_label = ref.get("region_label", "V3V4")
silva_version = str(ref.get("silva_version", "138.2"))
silva_target = ref.get("silva_target", "SSURef_NR99")
prim_fwd = ref.get("primers", {}).get("forward", "CCTACGGGNGGCWGCAG")
prim_rev = ref.get("primers", {}).get("reverse", "GACTACHVGGGTATCTAATCC")
symbiont_list = ref.get("symbiont_list", "refs/symbiont_taxa.txt")

SILVA_DIR = Path(f"refs/{region_label}-silva-{silva_version}")
SILVA_CLASSIFIER = SILVA_DIR / f"silva-{silva_version}-{region_label}-uniq-classifier.qza"
SILVA_UNIQ_SEQS = SILVA_DIR / f"silva-{silva_version}-{region_label}-uniq-seqs.qza"
SILVA_UNIQ_TAX  = SILVA_DIR / f"silva-{silva_version}-{region_label}-uniq-tax.qza"

gg2 = cfg.get("gg2", {})
GG2_SEQS = gg2.get("seqs_qza", "")
GG2_TAX  = gg2.get("tax_qza", "")
GG2_DIR  = Path(gg2.get("outdir", f"refs/GG2-{region_label}"))
GG2_CLASSIFIER = GG2_DIR / f"{region_label}-uniq-classifier.qza"

primary_classifier = Path(cfg.get("classifier_qza", str(SILVA_CLASSIFIER)))

conc = cfg.get("concordance", {})
TOP_N = int(conc.get("top_n", 50))

# ---- New: rank utilities for summaries ----
RANKS = {"phylum": 2, "family": 5, "genus": 6, "species": 7}
RANK_ORDER = ["kingdom","phylum","class","order","family","genus","species"]

rule all:
    input:
        expand("qiime2/qc/fastqc/{s}_R1_fastqc.html", s=SAMPLES),
        expand("qiime2/qc/fastqc/{s}_R2_fastqc.html", s=SAMPLES),
        "qiime2/qc/multiqc/multiqc_report.html",
        "qiime2/demux.qza",
        "qiime2/demux.qzv",
        "qiime2/table.qza",
        "qiime2/rep-seqs.qza",
        "qiime2/denoise-stats.qza",
        "qiime2/taxonomy/taxonomy_primary.qza",
        "qiime2/taxonomy/taxa-barplot_primary.qzv",
        expand("qiime2/taxonomy/taxonomy_gg2.qza", allow_missing=True),
        expand(f"qiime2/concordance/top{TOP_N}_concordance.csv", allow_missing=True)
        # 注：新加的“汇总表”不自动加进 all，以免改变你既有调用方式；按需单独触发即可。

# ---------------- QC ----------------
rule fastqc:
    input:
        r1 = lambda wc: RAW / f"{wc.sample}{R1SFX}",
        r2 = lambda wc: RAW / f"{wc.sample}{R2SFX}",
    output:
        html1 = "qiime2/qc/fastqc/{sample}_R1_fastqc.html",
        html2 = "qiime2/qc/fastqc/{sample}_R2_fastqc.html",
        zip1  = "qiime2/qc/fastqc/{sample}_R1_fastqc.zip",
        zip2  = "qiime2/qc/fastqc/{sample}_R2_fastqc.zip",
    conda: "workflow/envs/qc.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p qiime2/qc/fastqc
        fastqc -t {threads} -o qiime2/qc/fastqc {input.r1} {input.r2}
        """

rule multiqc:
    input:
        expand("qiime2/qc/fastqc/{s}_R1_fastqc.zip", s=SAMPLES),
        expand("qiime2/qc/fastqc/{s}_R2_fastqc.zip", s=SAMPLES),
    output:
        html = "qiime2/qc/multiqc/multiqc_report.html"
    conda: "workflow/envs/qc.yaml"
    shell:
        r"""
        mkdir -p qiime2/qc/multiqc
        multiqc -o qiime2/qc/multiqc qiime2/qc/fastqc
        """

# ------------- Cutadapt -------------
rule cutadapt_paired:
    input:
        r1 = lambda wc: RAW / f"{wc.sample}{R1SFX}",
        r2 = lambda wc: RAW / f"{wc.sample}{R2SFX}",
    output:
        r1t = "qiime2/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2t = "qiime2/trimmed/{sample}_R2.trimmed.fastq.gz"
    params:
        fwd = prim_fwd,
        rev = prim_rev,
        extra = cfg.get("cutadapt_extra","--minimum-length 50 -q 20,20")
    conda: "workflow/envs/cutadapt.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p qiime2/trimmed
        cutadapt -j {threads}           -g {params.fwd} -G {params.rev}           -o {output.r1t} -p {output.r2t}           {params.extra}           {input.r1} {input.r2} > qiime2/trimmed/{wildcards.sample}.cutadapt.log
        """

# ------------- Manifest -------------
rule make_manifest:
    input:
        r1s = expand("qiime2/trimmed/{s}_R1.trimmed.fastq.gz", s=SAMPLES),
        r2s = expand("qiime2/trimmed/{s}_R2.trimmed.fastq.gz", s=SAMPLES),
    output:
        manifest = "qiime2/manifest_pe.tsv"
    run:
        from pathlib import Path
        man = Path(output.manifest); man.parent.mkdir(parents=True, exist_ok=True)
        with man.open("w") as fh:
            fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
            for s in SAMPLES:
                f = (Path("qiime2/trimmed")/f"{s}_R1.trimmed.fastq.gz").resolve()
                r = (Path("qiime2/trimmed")/f"{s}_R2.trimmed.fastq.gz").resolve()
                fh.write(f"{s}\t{f}\t{r}\n")

# ------------- QIIME 2 core pipeline -------------
rule qiime_import:
    input: "qiime2/manifest_pe.tsv"
    output: "qiime2/demux.qza"
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    shell:
        r"""
        {Q2} tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path {input} --output-path {output}
        """

rule demux_summarize:
    input: "qiime2/demux.qza"
    output: "qiime2/demux.qzv"
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    shell:
        r"""
        {Q2} demux summarize --i-data {input} --o-visualization {output}
        """

rule dada2_denoise_paired:
    input: demux="qiime2/demux.qza"
    output:
        table="qiime2/table.qza",
        reps ="qiime2/rep-seqs.qza",
        stats="qiime2/denoise-stats.qza"
    params:
        tlf = cfg["dada2"].get("trim_left_f", 0),
        tlr = cfg["dada2"].get("trim_left_r", 0),
        trf = cfg["dada2"].get("trunc_len_f", 0),
        trr = cfg["dada2"].get("trunc_len_r", 0)
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    threads: 8
    shell:
        r"""
        {Q2} dada2 denoise-paired \
          --i-demultiplexed-seqs {input.demux} \
          --p-trim-left-f {params.tlf} --p-trim-left-r {params.tlr} \
          --p-trunc-len-f {params.trf} --p-trunc-len-r {params.trr} \
          --o-table {output.table} \
          --o-representative-sequences {output.reps} \
          --o-denoising-stats {output.stats} \
          --p-n-threads {threads}
        """

# ------------- Classifier build (SILVA) -------------
rule build_classifier_silva:
    output:
        classifier = str(SILVA_CLASSIFIER),
        uniq_seqs  = str(SILVA_UNIQ_SEQS),
        uniq_tax   = str(SILVA_UNIQ_TAX)
    params:
        region = region_label,
        ver = silva_version,
        target = silva_target,
        fwd = prim_fwd,
        rev = prim_rev,
        symb = symbiont_list
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    shell:
        r"""
        TARGET="{params.target}" \
        bash refs/get_silva_classifier.sh {params.region} '{params.fwd}' '{params.rev}' {params.ver} {params.symb}
        """

# ------------- Classifier build (GG2, 2024.09; optional) -------------
rule build_classifier_gg2:
    input:
        seqs = GG2_SEQS,
        tax  = GG2_TAX
    output:
        classifier = str(GG2_CLASSIFIER)
    params:
        region = region_label,
        outdir = str(GG2_DIR),
        fwd = prim_fwd,
        rev = prim_rev
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"

        # 1) extract V3–V4 reads from GG2 full-length backbone
        {Q2} feature-classifier extract-reads \
          --i-sequences "{input.seqs}" \
          --p-f-primer '{params.fwd}' \
          --p-r-primer '{params.rev}' \
          --o-reads "{params.outdir}/reads.qza"

        # 2) RESCRIPt dereplicate (auto-pick taxonomy flag name)
        DEREP_FLAG="--o-dereplicated-taxa"
        {Q2} rescript dereplicate --help | grep -q -- "--o-dereplicated-taxonomy" && DEREP_FLAG="--o-dereplicated-taxonomy"

        {Q2} rescript dereplicate \
          --i-sequences "{params.outdir}/reads.qza" \
          --i-taxa      "{input.tax}" \
          --p-mode uniq \
          --o-dereplicated-sequences "{params.outdir}/uniq-seqs.qza" \
          $DEREP_FLAG                     "{params.outdir}/uniq-tax.qza"

        # 3) train Naive Bayes classifier
        {Q2} feature-classifier fit-classifier-naive-bayes \
          --i-reference-reads    "{params.outdir}/uniq-seqs.qza" \
          --i-reference-taxonomy "{params.outdir}/uniq-tax.qza" \
          --o-classifier         "{output.classifier}"
        """

# ------------- Taxonomy (Primary) -------------
rule taxonomy_primary:
    input:
        reps="qiime2/rep-seqs.qza",
        classifier=str(primary_classifier),
        meta=str(META)
    output:
        tax="qiime2/taxonomy/taxonomy_primary.qza",
        bar="qiime2/taxonomy/taxa-barplot_primary.qzv"
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    shell:
        r"""
        mkdir -p qiime2/taxonomy
        {Q2} feature-classifier classify-sklearn --i-classifier {input.classifier} --i-reads {input.reps} --o-classification {output.tax}
        {Q2} taxa barplot --i-table qiime2/table.qza --i-taxonomy {output.tax} --m-metadata-file "{input.meta}" --o-visualization {output.bar}
        """

# ------------- Taxonomy (GG2 cross-check, optional) -------------
use_gg2 = bool(GG2_SEQS) and bool(GG2_TAX)
rule export_taxonomy_primary_tsv:
    input: "qiime2/taxonomy/taxonomy_primary.qza"
    output: tsv="qiime2/taxonomy/taxonomy_primary.tsv"
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
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
    conda: None
    shell:
        r"""
        rm -rf qiime2/_tmp || true; mkdir -p qiime2/_tmp
        {Q2} tools export --input-path {input} --output-path qiime2/_tmp
        conda run -n {QIIME_ENV_NAME} biom convert -i {output.biom} -o {output.tsv} --to-tsv
        """	    	
		
if use_gg2:
    rule taxonomy_gg2:
        input:
            reps="qiime2/rep-seqs.qza",
            classifier=str(GG2_CLASSIFIER)
        output:
            tax="qiime2/taxonomy/taxonomy_gg2.qza"
        conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
        shell:
            r"""
            mkdir -p qiime2/taxonomy
            {Q2} feature-classifier classify-sklearn --i-classifier {input.classifier} --i-reads {input.reps} --o-classification {output.tax}
            """
   
    rule export_taxonomy_gg2_tsv:
        input: "qiime2/taxonomy/taxonomy_gg2.qza"
        output: tsv="qiime2/taxonomy/taxonomy_gg2.tsv"
        conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
        shell:
            r"""
            rm -rf qiime2/taxonomy/_export_gg2 || true
            {Q2} tools export --input-path {input} --output-path qiime2/taxonomy/_export_gg2
            mv qiime2/taxonomy/_export_gg2/taxonomy.tsv {output.tsv}
            rm -rf qiime2/taxonomy/_export_gg2
            """

    # --- ASV (GG2 taxonomy) → 宽表 ---
    rule asv_level_table_gg2:
        input:
            table_qza="qiime2/table.qza",
            tax_qza="qiime2/taxonomy/taxonomy_gg2.qza"
        output:
            tsv = "qiime2/summary/asv_gg2_{rank,phylum|family|genus|species}.tsv"
        params:
            rank = lambda wc: wc.rank  # phylum/family/genus/species
        conda: None
        shell:
            r"""
            set -euo pipefail
            TMP=qiime2/summary/_tmp_asv_gg2_{wildcards.rank}
            mkdir -p "$TMP/tab" "$TMP/tax"
            {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
            conda run -n {QIIME_ENV_NAME} biom convert \
              -i "$TMP/tab/feature-table.biom" \
              -o "$TMP/feature-table.tsv" --to-tsv
            {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

            conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
              --table-tsv "$TMP/feature-table.tsv" \
              --tax-tsv   "$TMP/tax/taxonomy.tsv" \
              --rank {params.rank} \
              --output {output.tsv}
            """

    # --- ASV (GG2) → 相对丰度 + 长表 ---
    rule asv_gg2_post_rank:
        input:
            wide = "qiime2/summary/asv_gg2_{rank}.tsv",
            meta = str(META)
        output:
            rel  = "qiime2/summary/asv_gg2_{rank}.relabund.tsv",
            long = "qiime2/summary/asv_gg2_{rank}.long.tsv"
        params:
            group = cfg.get("summary", {}).get("group_column", "group")
        conda: None
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
              --wide-tsv {input.wide} \
              --metadata {input.meta} \
              --group-col {params.group} \
              --out-rel  {output.rel} \
              --out-long {output.long}
            """
            
    # OTU97 代表序列用 GG2 分类
    rule otu97_taxonomy_gg2:
        input:
            reps="qiime2/otu97/otu97-rep-seqs.qza",
            classifier=str(GG2_CLASSIFIER)   # 已在前面定义
        output:
            tax="qiime2/otu97/taxonomy_otu97_gg2.qza"
        conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
        shell:
            r"""
            {Q2} feature-classifier classify-sklearn \
              --i-classifier {input.classifier} \
              --i-reads {input.reps} \
              --o-classification {output.tax}
            """

    # OTU97 (GG2) → 宽表
    rule otu97_gg2_level_table:
        input:
            table_qza="qiime2/otu97/otu97-table.qza",
            tax_qza  ="qiime2/otu97/taxonomy_otu97_gg2.qza"
        output:
            tsv="qiime2/summary/otu97_gg2_{rank,phylum|family|genus|species}.tsv"
        params:
            rank=lambda wc: wc.rank
        conda: None
        shell:
            r"""
            set -euo pipefail
            TMP=qiime2/summary/_tmp_otu97_gg2_{wildcards.rank}
            mkdir -p "$TMP/tab" "$TMP/tax"
            {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
            conda run -n {QIIME_ENV_NAME} biom convert \
              -i "$TMP/tab/feature-table.biom" \
              -o "$TMP/otu-table.tsv" --to-tsv
            {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

            conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
              --table-tsv "$TMP/otu-table.tsv" \
              --tax-tsv   "$TMP/tax/taxonomy.tsv" \
              --rank {params.rank} \
              --output {output.tsv}
            """

    # OTU97 (GG2) → 相对丰度 + 长表
    rule otu97_gg2_post_rank:
        input:
            wide = "qiime2/summary/otu97_gg2_{rank}.tsv",
            meta = str(META)
        output:
            rel  = "qiime2/summary/otu97_gg2_{rank}.relabund.tsv",
            long = "qiime2/summary/otu97_gg2_{rank}.long.tsv"
        params:
            group = cfg.get("summary", {}).get("group_column", "group")
        conda: None
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
              --wide-tsv {input.wide} \
              --metadata {input.meta} \
              --group-col {params.group} \
              --out-rel  {output.rel} \
              --out-long {output.long}
            """
 
    # 1) 并列注释清单 + 各阶元一致性统计
    rule pair_taxonomy_tables:
        input:
            silva="qiime2/taxonomy/taxonomy_primary.tsv",
            gg2  ="qiime2/taxonomy/taxonomy_gg2.tsv"
        output:
            pairs ="qiime2/concordance/pair_taxonomy.tsv",
            agree ="qiime2/concordance/agreement_by_rank.tsv"
        conda: None
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/pair_tax.py \
              --silva {input.silva} \
              --gg2   {input.gg2} \
              --out-pairs {output.pairs} \
              --out-agree {output.agree}
            """

    # 2) （可选）按“最深一致阶元”产生 resolved 注释（QIIME TSV 格式）
    rule resolved_taxonomy_tsv:
        input:
            silva="qiime2/taxonomy/taxonomy_primary.tsv",
            gg2  ="qiime2/taxonomy/taxonomy_gg2.tsv"
        output:
            tsv="qiime2/taxonomy/taxonomy_resolved.tsv"
        conda: None
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/resolve_tax.py \
              --silva {input.silva} \
              --gg2   {input.gg2} \
              --out   {output.tsv}
            """

    # 3) 用 resolved 注释生成各阶元宽表（复用 rank_wide.py；从 TSV 直接读注释）
    wildcard_constraints:
        rrank="phylum|family|genus|species"

    rule asv_level_table_resolved:
        input:
            table_qza="qiime2/table.qza",
            tax_tsv  ="qiime2/taxonomy/taxonomy_resolved.tsv"
        output:
            tsv="qiime2/summary/asv_resolved_{rrank}.tsv"
        params:
            rrank=lambda wc: wc.rrank
        conda: None
        shell:
            r"""
            set -euo pipefail
            TMP=qiime2/summary/_tmp_asv_resolved_{wildcards.rrank}
            mkdir -p "$TMP/tab"
            {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
            conda run -n {QIIME_ENV_NAME} biom convert \
              -i "$TMP/tab/feature-table.biom" \
              -o "$TMP/feature-table.tsv" --to-tsv

            conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
              --table-tsv "$TMP/feature-table.tsv" \
              --tax-tsv   {input.tax_tsv} \
              --rank      {params.rrank} \
              --output    {output.tsv}
            """

    # （可选）resolved 的相对丰度/长表
    rule asv_resolved_post_rank:
        input:
            wide = "qiime2/summary/asv_resolved_{rrank}.tsv",
            meta = str(META)
        output:
            rel  = "qiime2/summary/asv_resolved_{rrank}.relabund.tsv",
            long = "qiime2/summary/asv_resolved_{rrank}.long.tsv"
        params:
            group = cfg.get("summary", {}).get("group_column", "group")
        conda: None
        shell:
            r"""
            conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
              --wide-tsv {input.wide} \
              --metadata {input.meta} \
              --group-col {params.group} \
              --out-rel  {output.rel} \
              --out-long {output.long}
            """
 
    rule concordance_top:
        input:
            table_tsv="qiime2/feature-table.tsv",
            silva_tsv="qiime2/taxonomy/taxonomy_primary.tsv",
            gg2_tsv="qiime2/taxonomy/taxonomy_gg2.tsv"
        output:
            csv=f"qiime2/concordance/top{TOP_N}_concordance.csv"
        conda: "workflow/envs/qc.yaml"
        script:
            "workflow/scripts/concordance.py"

# ---------------- 97% OTU (de novo clustering) ----------------
rule otu97_denovo:
    input:
        table="qiime2/table.qza",
        reps ="qiime2/rep-seqs.qza"
    output:
        t="qiime2/otu97/otu97-table.qza",
        r="qiime2/otu97/otu97-rep-seqs.qza"
    threads: 8
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    shell:
        r"""
        mkdir -p qiime2/otu97
        {Q2} vsearch cluster-features-de-novo \
          --i-table {input.table} \
          --i-sequences {input.reps} \
          --p-perc-identity 0.97 \
          --o-clustered-table {output.t} \
          --o-clustered-sequences {output.r}
        """

# -------- OTU97 taxonomy with primary classifier (SILVA) --------
rule otu97_taxonomy_primary:
    input:
        reps="qiime2/otu97/otu97-rep-seqs.qza",
        classifier=str(primary_classifier)
    output:
        tax="qiime2/otu97/taxonomy_otu97_primary.qza"
    conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
    shell:
        r"""
        {Q2} feature-classifier classify-sklearn \
          --i-classifier {input.classifier} \
          --i-reads {input.reps} \
          --o-classification {output.tax}
        """

# ------------- Phylogeny via SEPP (GG2 non-V4 backbone) -------------
sepp_ok = bool(GG2_SEPP_REF)

if sepp_ok:
    rule gg2_nonv4_sepp_insertion:
        input:
            reps = "qiime2/rep-seqs.qza",
            ref  = str(GG2_SEPP_REF)
        output:
            tree = "qiime2/trees/gg2_nonv4_sepp/insertion-tree.qza",
            plc  = "qiime2/trees/gg2_nonv4_sepp/placements.qza"
        threads: 8
        conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
        shell:
            r"""
            set -euo pipefail
            mkdir -p qiime2/trees/gg2_nonv4_sepp
            {Q2} fragment-insertion sepp \
              --i-representative-sequences {input.reps} \
              --i-reference-database      {input.ref}  \
              --p-threads {threads} \
              --o-tree       {output.tree} \
              --o-placements {output.plc}
            """

    rule export_sepp_tree_newick:
        input: "qiime2/trees/gg2_nonv4_sepp/insertion-tree.qza"
        output: "qiime2/trees/gg2_nonv4_sepp/insertion-tree.nwk"
        conda: "workflow/envs/qiime2.yaml" if not USE_EXISTING_QIIME else None
        shell:
            r"""
            {Q2} tools export --input-path {input} --output-path qiime2/trees/gg2_nonv4_sepp/_export
            mv qiime2/trees/gg2_nonv4_sepp/_export/tree.nwk {output}
            rm -rf qiime2/trees/gg2_nonv4_sepp/_export
            """

# =====================================================================
#                 ASV & OTU SUMMARY WIDE TABLES
# =====================================================================
# --- ASV summaries per rank (phylum/family/genus/species) ---
# 直接从 table.qza + taxonomy_primary.qza 自导出临时 TSV，再汇总成宽表。
rule asv_level_table:
    input:
        table_qza="qiime2/table.qza",
        tax_qza="qiime2/taxonomy/taxonomy_primary.qza"
    output:
        tsv = "qiime2/summary/asv_{rank,phylum|family|genus|species}.tsv"
    params:
        rank = lambda wc: wc.rank  # phylum/family/genus/species
    conda: None
    shell:
        r"""
        set -euo pipefail
        TMP=qiime2/summary/_tmp_asv_{wildcards.rank}
        mkdir -p "$TMP/tab" "$TMP/tax"
        # 导出 QZA → TSV
        {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
        conda run -n {QIIME_ENV_NAME} biom convert \
          -i "$TMP/tab/feature-table.biom" \
          -o "$TMP/feature-table.tsv" --to-tsv
        {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

        # 用独立脚本生成宽表
        conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
          --table-tsv "$TMP/feature-table.tsv" \
          --tax-tsv "$TMP/tax/taxonomy.tsv" \
          --rank {params.rank} \
          --output {output.tsv}
        """

# --- 97% OTU (de novo) + taxonomy + per-rank summaries ---
rule otu97_level_table:
    input:
        table_qza="qiime2/otu97/otu97-table.qza",
        tax_qza  ="qiime2/otu97/taxonomy_otu97_primary.qza"
    output:
        tsv="qiime2/summary/otu97_{rank,phylum|family|genus|species}.tsv"
    params:
        rank=lambda wc: wc.rank
    conda: None
    shell:
        r"""
        set -euo pipefail
        TMP=qiime2/summary/_tmp_otu97_{wildcards.rank}
        mkdir -p "$TMP/tab" "$TMP/tax"
        {Q2} tools export --input-path {input.table_qza} --output-path "$TMP/tab"
        conda run -n {QIIME_ENV_NAME} biom convert \
          -i "$TMP/tab/feature-table.biom" \
          -o "$TMP/otu-table.tsv" --to-tsv
        {Q2} tools export --input-path {input.tax_qza} --output-path "$TMP/tax"

        conda run -n {QIIME_ENV_NAME} python workflow/scripts/rank_wide.py \
          --table-tsv "$TMP/otu-table.tsv" \
          --tax-tsv "$TMP/tax/taxonomy.tsv" \
          --rank {params.rank} \
          --output {output.tsv}
        """

# --- ASV → 相对丰度 + 长表 ---
rule asv_post_rank:
    input:
        wide = "qiime2/summary/asv_{rank}.tsv",
        meta = str(META)
    output:
        rel  = "qiime2/summary/asv_{rank}.relabund.tsv",
        long = "qiime2/summary/asv_{rank}.long.tsv"
    params:
        group = cfg.get("summary", {}).get("group_column", "group")  # 按需在 config.yaml 里设定
    conda: None
    shell:
        r"""
        conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
          --wide-tsv {input.wide} \
          --metadata {input.meta} \
          --group-col {params.group} \
          --out-rel  {output.rel} \
          --out-long {output.long}
        """

# --- OTU97 → 相对丰度 + 长表 ---
rule otu97_post_rank:
    input:
        wide = "qiime2/summary/otu97_{rank}.tsv",
        meta = str(META)
    output:
        rel  = "qiime2/summary/otu97_{rank}.relabund.tsv",
        long = "qiime2/summary/otu97_{rank}.long.tsv"
    params:
        group = cfg.get("summary", {}).get("group_column", "group")
    conda: None
    shell:
        r"""
        conda run -n {QIIME_ENV_NAME} python workflow/scripts/wide_post.py \
          --wide-tsv {input.wide} \
          --metadata {input.meta} \
          --group-col {params.group} \
          --out-rel  {output.rel} \
          --out-long {output.long}
        """
