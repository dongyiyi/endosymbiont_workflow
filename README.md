# endosymbiont workflow [![CI](https://github.com/dongyiyi/endosymbiont_workflow/actions/workflows/ci.yml/badge.svg)](https://github.com/dongyiyi/endosymbiont_workflow/actions/workflows/ci.yml)

Endosymbiont 16S V3–V4 Workflow (**Snakemake + QIIME 2 [https://qiime2.org/]**)
Reproducible pipeline for bark & ambrosia beetle **endosymbiont** 16S amplicons.
Raw FASTQs → QC/denoise (DADA2) → taxonomy (**SILVA [https://www.arb-silva.de/]** primary; optional **Greengenes2（GG2）[https://greengenes2.ucsd.edu/]** cross-check) → ASV/OTU97 wide/long/rel-abundance tables → optional phylogeny.


## 1) Requirements

**Conda/Mamba** (recommended)

**Snakemake 8** (driver env)

**QIIME 2 amplicon 2024.10** (analysis env)

Create driver env:
```
conda create -n snakemake -c conda-forge -c bioconda snakemake=8 python=3.11 graphviz -y
```

Use your normal method to install/activate **qiime2-amplicon-2024.10** (the workflow assumes this env name; changeable via `workflow/config.yaml`).

## 2) Project layout
```
project/
├─ Snakefile
├─ workflow/
│  ├─ config.yaml
│  ├─ envs/                  # conda env specs
│  └─ scripts/               # rank_wide.py, wide_post.py, (optional) resolve_tax.py
├─ raw/                      # SampleX_R1.fastq.gz, SampleX_R2.fastq.gz
├─ metadata/
│  └─ sample-metadata.tsv    # QIIME 2 metadata (TSV; first col: sample-id)
└─ refs/                     # reference DBs, trained classifiers, helper lists
```

**Metadata gotcha**: the first header must be `sample-id` (exact), file must be **tab-separated**, and sample IDs must match FASTQ prefixes.

## 3) Configure `workflow/config.yaml`
```
raw_dir: raw
metadata_tsv: metadata/sample-metadata.tsv

r1_suffix: _R1.fastq.gz
r2_suffix: _R2.fastq.gz

cutadapt_extra: --minimum-length 50 -q 20,20

dada2:
  trim_left_f: 0
  trim_left_r: 0
  trunc_len_f: 0
  trunc_len_r: 0

use_existing_qiime_env: true
qiime_env_name: qiime2-amplicon-2024.10

### primary (SILVA) classifier for V3–V4
classifier_qza: refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-classifier.qza

reference:
  region_label: V3V4
  silva_version: '138.2'
  silva_target: SSURef_NR99
  primers:
    forward: CCTACGGGNGGCWGCAG
    reverse: GACTACHVGGGTATCTAATCC
  symbiont_list: refs/symbiont_taxa.txt

gg2:
  seqs_qza: refs/GG2-raw/2024.09.backbone.full-length.fna.qza
  tax_qza:  refs/GG2-raw/2024.09.backbone.tax.qza
  outdir:   refs/GG2-V3V4
  sepp_ref_qza: refs/GG2-raw/2022.10.backbone.sepp-reference.qza

summary:
  group_column: group         # e.g. group, native_status, region
concordance:
  top_n: 50
```
## 4) Reference databases

**Train once per primer set.** If future samples use the same V3–V4 primers, reuse the trained classifiers—no need to rebuild.

### 4.1 SILVA 138.2 (SSURef NR99) → V3–V4 classifier

Start from ready QZA (recommended) or import from FASTA/TSV.

Import (if needed)
```
# sequences
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path silva-138.2-ssu-nr99-seqs.fasta \
  --output-path refs/V3V4-silva-138.2/silva-138.2-ssu-nr99-seqs.qza

# taxonomy (2 cols: FeatureID<TAB>Taxon)
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path silva-138.2-ssu-nr99-tax.tsv \
  --output-path refs/V3V4-silva-138.2/silva-138.2-ssu-nr99-tax.qza
```

**Extract V3–V4, dereplicate, train**
```
# 1) extract reads with your primers
qiime rescript extract-reads \
  --i-sequences refs/V3V4-silva-138.2/silva-138.2-ssu-nr99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-n-jobs 8 \
  --o-reads   refs/V3V4-silva-138.2/silva-138.2-V3V4-reads.qza

# 2) dereplicate with taxonomy
qiime rescript dereplicate \
  --i-sequences refs/V3V4-silva-138.2/silva-138.2-V3V4-reads.qza \
  --i-taxa      refs/V3V4-silva-138.2/silva-138.2-ssu-nr99-tax.qza \
  --p-mode uniq \
  --o-dereplicated-sequences refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-seqs.qza \
  --o-dereplicated-taxa      refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-tax.qza

# 3) train classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads     refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-seqs.qza \
  --i-reference-taxonomy  refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-tax.qza \
  --o-classifier          refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-classifier.qza

qiime tools validate refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-classifier.qza
```
### 4.2 GG2 backbone → V3–V4 classifier

**Import (if needed)**
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 2024.09.backbone.full-length.fna \
  --output-path refs/GG2-raw/2024.09.backbone.full-length.fna.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 2024.09.backbone.tax.tsv \
  --output-path refs/GG2-raw/2024.09.backbone.tax.qza
```

**Extract, dereplicate, train**
```
mkdir -p refs/GG2-V3V4

qiime rescript extract-reads \
  --i-sequences refs/GG2-raw/2024.09.backbone.full-length.fna.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-n-jobs 8 \
  --o-reads   refs/GG2-V3V4/reads.qza

qiime rescript dereplicate \
  --i-sequences refs/GG2-V3V4/reads.qza \
  --i-taxa      refs/GG2-raw/2024.09.backbone.tax.qza \
  --p-mode uniq \
  --o-dereplicated-sequences refs/GG2-V3V4/V3V4-uniq-seqs.qza \
  --o-dereplicated-taxa      refs/GG2-V3V4/V3V4-uniq-tax.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads     refs/GG2-V3V4/V3V4-uniq-seqs.qza \
  --i-reference-taxonomy  refs/GG2-V3V4/V3V4-uniq-tax.qza \
  --o-classifier          refs/GG2-V3V4/V3V4-uniq-classifier.qza

qiime tools validate refs/GG2-V3V4/V3V4-uniq-classifier.qza
```
### 4.3 Why `symbiont_taxa.txt`?

A curated list of known **endosymbiont** clades (e.g., genera/families reported in beetles). We include it to:

1. Optionally build a **symbiont-focused classifier** (subset of SILVA) to increase precision/sensitivity for target lineages.

2. Quickly **subset downstream** tables to symbiont rows for focused visualization/statistics.

(We still run the **full** SILVA classifier as the primary taxonomy; the focused version is an additional check.)

**Build focused classifier (optional)**
```
INCLUDE=$(paste -sd, refs/symbiont_taxa.txt)

qiime taxa filter-seqs \
  --i-sequences refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-seqs.qza \
  --i-taxonomy refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-tax.qza \
  --p-include "$INCLUDE" \
  --o-filtered-sequences refs/V3V4-silva-138.2/symbiont-V3V4-uniq-seqs.qza

qiime taxa filter-taxonomy \
  --i-taxonomy refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-tax.qza \
  --p-include "$INCLUDE" \
  --o-filtered-taxonomy refs/V3V4-silva-138.2/symbiont-V3V4-uniq-tax.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads     refs/V3V4-silva-138.2/symbiont-V3V4-uniq-seqs.qza \
  --i-reference-taxonomy  refs/V3V4-silva-138.2/symbiont-V3V4-uniq-tax.qza \
  --o-classifier          refs/V3V4-silva-138.2/symbiont-V3V4-uniq-classifier.qza
```
## 5) Run the workflow

Activate `snakemake` env when running commands below.

### 5.1 Dry-run / summary
```
conda run -n snakemake snakemake \
  -npr --summary \
  --configfile workflow/config.yaml
```
### 5.2 Core analysis (SILVA)
```
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/demux.qzv \
  qiime2/table.qza \
  qiime2/rep-seqs.qza \
  qiime2/taxonomy/taxa-barplot_primary.qzv
```

Open `.qzv` locally with **QIIME 2 View** (or download to your workstation and drag-drop).

### 5.3 ASV/OTU wide/long/rel-abundance tables (SILVA)
```
# ASV (genus & species, wide + long + relabund)
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/summary/asv_genus.tsv \
  qiime2/summary/asv_genus.long.tsv \
  qiime2/summary/asv_genus.relabund.tsv \
  qiime2/summary/asv_species.tsv \
  qiime2/summary/asv_species.long.tsv \
  qiime2/summary/asv_species.relabund.tsv

# OTU97 and their summaries
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/otu97/otu97-table.qza \
  qiime2/otu97/otu97-rep-seqs.qza \
  qiime2/summary/otu97_genus.tsv \
  qiime2/summary/otu97_genus.long.tsv \
  qiime2/summary/otu97_genus.relabund.tsv \
  qiime2/summary/otu97_species.tsv \
  qiime2/summary/otu97_species.long.tsv \
  qiime2/summary/otu97_species.relabund.tsv
```
### 5.4 GG2 cross-check (optional)
```
# taxonomy and exports
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/taxonomy/taxonomy_gg2.qza \
  qiime2/taxonomy/taxonomy_gg2.tsv

# ASV wide/long/relabund using GG2 taxonomy
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/summary/asv_gg2_genus.tsv \
  qiime2/summary/asv_gg2_genus.long.tsv \
  qiime2/summary/asv_gg2_genus.relabund.tsv \
  qiime2/summary/asv_gg2_species.tsv \
  qiime2/summary/asv_gg2_species.long.tsv \
  qiime2/summary/asv_gg2_species.relabund.tsv

# OTU97 GG2 summaries
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/otu97/taxonomy_otu97_gg2.qza \
  qiime2/summary/otu97_gg2_genus.tsv \
  qiime2/summary/otu97_gg2_genus.long.tsv \
  qiime2/summary/otu97_gg2_genus.relabund.tsv \
  qiime2/summary/otu97_gg2_species.tsv \
  qiime2/summary/otu97_gg2_species.long.tsv \
  qiime2/summary/otu97_gg2_species.relabund.tsv
```
### 5.5 SILVA vs GG2 concordance
```
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/concordance/top50_concordance.csv
```
### 5.6 Optional phylogeny (SEPP placement)

Requires `gg2.sepp_ref_qza` in config:
```
conda run -n snakemake snakemake \
  --configfile workflow/config.yaml --cores 8 --use-conda --printshellcmds \
  qiime2/trees/gg2_nonv4_sepp/insertion-tree.qza \
  qiime2/trees/gg2_nonv4_sepp/placements.qza \
  qiime2/trees/gg2_nonv4_sepp/insertion-tree.nwk
```

## 6) One-command “run all” script

Create `run_all.sh` and make it executable (`chmod +x run_all.sh`):
```
#!/usr/bin/env bash
#SBATCH --job-name=endsymbiont_runall
#SBATCH --output=logs/run_all.%j.out
#SBATCH --error=logs/run_all.%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# -----------------------------
# Config (edit if needed)
# -----------------------------
SNAKEMAKE_ENV="snakemake"
QIIME_ENV="qiime2-amplicon-2024.10"
R_ENV="r-microbiome"

CONFIG="workflow/config.yaml"
SNAKE1="Snakefile"
SNAKE2="Snakefile_v2"

META_BASIC="metadata/sample-metadata.tsv"
META_DETAIL="metadata/sample-metadata_details.tsv"

# Sampling depth for rarefaction/core-metrics (your final choice)
SAMPLING_DEPTH="500000"

# Where to store outputs
mkdir -p logs

# Use Slurm cores if available, otherwise fallback
CORES="${SLURM_CPUS_PER_TASK:-8}"

echo "[INFO] Start run_all.sh at: $(date)"
echo "[INFO] Using CORES=${CORES}"
echo "[INFO] Using SAMPLING_DEPTH=${SAMPLING_DEPTH}"

# -----------------------------
# Helper: run snakemake safely
# -----------------------------
run_snakemake () {
  local snakefile="$1"
  echo "[INFO] Running Snakemake with ${snakefile}"
  conda run -n "${SNAKEMAKE_ENV}" snakemake \
    --snakefile "${snakefile}" \
    --configfile "${CONFIG}" \
    --cores "${CORES}" \
    --use-conda \
    --printshellcmds \
    --rerun-incomplete \
    --rerun-triggers mtime
}

# -----------------------------
# 1) Main pipeline (Snakefile)
# -----------------------------
if [[ -f "${SNAKE1}" ]]; then
  run_snakemake "${SNAKE1}"
else
  echo "[WARN] ${SNAKE1} not found; skipping."
fi

# -----------------------------
# 2) Additional/updated rules (Snakefile_v2)
# -----------------------------
if [[ -f "${SNAKE2}" ]]; then
  run_snakemake "${SNAKE2}"
else
  echo "[WARN] ${SNAKE2} not found; skipping."
fi

# -----------------------------
# 3) Post-QIIME2 summaries you ran manually
#    - These do NOT affect existing results unless output filenames collide.
# -----------------------------
mkdir -p qiime2/summary

# DADA2 stats tabulate/export (safe: writes new qzv folder)
if [[ -f "qiime2/denoise-stats.qza" ]]; then
  echo "[INFO] QIIME2: metadata tabulate denoise-stats"
  conda run -n "${QIIME_ENV}" qiime metadata tabulate \
    --m-input-file qiime2/denoise-stats.qza \
    --o-visualization qiime2/summary/denoise-stats.qzv

  echo "[INFO] QIIME2: export denoise-stats"
  rm -rf qiime2/summary/_export_denoise_stats || true
  conda run -n "${QIIME_ENV}" qiime tools export \
    --input-path qiime2/denoise-stats.qza \
    --output-path qiime2/summary/_export_denoise_stats
else
  echo "[WARN] qiime2/denoise-stats.qza not found; skip denoise-stats export."
fi

# Feature-table summary (safe: new qzv)
if [[ -f "qiime2/table.qza" ]]; then
  echo "[INFO] QIIME2: feature-table summarize"
  # Prefer detail metadata if exists; fallback to basic
  META_USE="${META_DETAIL}"
  if [[ ! -f "${META_USE}" ]]; then META_USE="${META_BASIC}"; fi

  conda run -n "${QIIME_ENV}" qiime feature-table summarize \
    --i-table qiime2/table.qza \
    --m-sample-metadata-file "${META_USE}" \
    --o-visualization qiime2/summary/table_summary.qzv
else
  echo "[WARN] qiime2/table.qza not found; skip feature-table summarize."
fi

# -----------------------------
# 4) Phylogeny (align-to-tree-mafft-fasttree)
#    - Avoid re-running if all outputs validate.
# -----------------------------
mkdir -p qiime2/trees

need_phylo="0"
for f in qiime2/trees/aligned-rep-seqs.qza \
         qiime2/trees/masked-aligned-rep-seqs.qza \
         qiime2/trees/unrooted-tree.qza \
         qiime2/trees/rooted-tree.qza
do
  if [[ ! -f "$f" ]]; then
    need_phylo="1"
  fi
done

if [[ "${need_phylo}" == "1" ]]; then
  if [[ -f "qiime2/rep-seqs.qza" ]]; then
    echo "[INFO] QIIME2: align-to-tree-mafft-fasttree (this can take long; run on compute node)"
    conda run -n "${QIIME_ENV}" qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences qiime2/rep-seqs.qza \
      --o-alignment qiime2/trees/aligned-rep-seqs.qza \
      --o-masked-alignment qiime2/trees/masked-aligned-rep-seqs.qza \
      --o-tree qiime2/trees/unrooted-tree.qza \
      --o-rooted-tree qiime2/trees/rooted-tree.qza \
      --p-n-threads "${CORES}"
  else
    echo "[WARN] qiime2/rep-seqs.qza not found; skip phylogeny build."
  fi
else
  echo "[INFO] Phylogeny outputs already exist; validating..."
  conda run -n "${QIIME_ENV}" qiime tools validate qiime2/trees/aligned-rep-seqs.qza || echo "[WARN] aligned invalid"
  conda run -n "${QIIME_ENV}" qiime tools validate qiime2/trees/masked-aligned-rep-seqs.qza || echo "[WARN] masked invalid"
  conda run -n "${QIIME_ENV}" qiime tools validate qiime2/trees/unrooted-tree.qza || echo "[WARN] unrooted invalid"
  conda run -n "${QIIME_ENV}" qiime tools validate qiime2/trees/rooted-tree.qza || echo "[WARN] rooted invalid"
fi

# -----------------------------
# 5) Alpha rarefaction (500k) + export HTML
# -----------------------------
if [[ -f "qiime2/table.qza" && -f "qiime2/trees/rooted-tree.qza" ]]; then
  META_USE="${META_DETAIL}"
  if [[ ! -f "${META_USE}" ]]; then META_USE="${META_BASIC}"; fi

  echo "[INFO] QIIME2: alpha-rarefaction at depth ${SAMPLING_DEPTH}"
  conda run -n "${QIIME_ENV}" qiime diversity alpha-rarefaction \
    --i-table qiime2/table.qza \
    --i-phylogeny qiime2/trees/rooted-tree.qza \
    --m-metadata-file "${META_USE}" \
    --p-max-depth "${SAMPLING_DEPTH}" \
    --o-visualization "qiime2/summary/alpha_rarefaction_${SAMPLING_DEPTH}.qzv"

  echo "[INFO] Export alpha-rarefaction qzv -> HTML"
  rm -rf "qiime2/summary/export_alpha_rarefaction_${SAMPLING_DEPTH}" || true
  conda run -n "${QIIME_ENV}" qiime tools export \
    --input-path "qiime2/summary/alpha_rarefaction_${SAMPLING_DEPTH}.qzv" \
    --output-path "qiime2/summary/export_alpha_rarefaction_${SAMPLING_DEPTH}"
else
  echo "[WARN] Missing table/rooted tree; skip alpha-rarefaction."
fi

# -----------------------------
# 6) Core metrics phylogenetic (500k)
# -----------------------------
if [[ -f "qiime2/table.qza" && -f "qiime2/trees/rooted-tree.qza" ]]; then
  META_USE="${META_DETAIL}"
  if [[ ! -f "${META_USE}" ]]; then META_USE="${META_BASIC}"; fi

  OUT_CORE="qiime2/diversity/core_phylo_${SAMPLING_DEPTH}"
  if [[ ! -d "${OUT_CORE}" ]]; then
    mkdir -p "$(dirname "${OUT_CORE}")"
    echo "[INFO] QIIME2: core-metrics-phylogenetic at depth ${SAMPLING_DEPTH}"
    conda run -n "${QIIME_ENV}" qiime diversity core-metrics-phylogenetic \
      --i-table qiime2/table.qza \
      --i-phylogeny qiime2/trees/rooted-tree.qza \
      --m-metadata-file "${META_USE}" \
      --p-sampling-depth "${SAMPLING_DEPTH}" \
      --output-dir "${OUT_CORE}"
  else
    echo "[INFO] ${OUT_CORE} already exists; skip core-metrics."
  fi
else
  echo "[WARN] Missing table/rooted tree; skip core-metrics-phylogenetic."
fi

# -----------------------------
# 7) R analyses (symbionts) you are running
#    - Run only if scripts exist.
# -----------------------------
run_R () {
  local script="$1"
  if [[ -f "${script}" ]]; then
    echo "[INFO] Running R script: ${script}"
    conda run -n "${R_ENV}" Rscript "${script}"
  else
    echo "[WARN] R script not found: ${script}; skipping."
  fi
}

run_R "qiime2/analysis_symbiont/run_symbiont_region_ploidy.R"
run_R "qiime2/analysis_symbiont/run_symbiont_presence_absence_extra.R"
run_R "qiime2/analysis_symbiont/run_symbiont_presence_absence_models_v2.R"

# (Optional) if you created these additional scripts
run_R "qiime2/analysis_symbiont/run_symbiont_species_level.R"
run_R "qiime2/analysis_symbiont/run_symbiont_conspecific_pairs.R"
run_R "qiime2/analysis_symbiont/run_symbiont_presence_threshold_sensitivity.R"

# -----------------------------
# 8) Package results for download
#    - Keep original folder
#    - Exclude fastq and fastqc.zip as requested
# -----------------------------
# echo "[INFO] Packaging qiime2 -> qiime2_results.gz (excluding fastq.gz and fastqc.zip)"
# tar -czf qiime2_results.gz \
#   --exclude='*.fastq.gz' \
#   --exclude='*fastqc.zip' \
#   qiime2

# Optional extra excludes (uncomment if you want smaller archive)
#   --exclude='qiime2/demux*' \
#   --exclude='qiime2/*/cutadapt*' \
#   --exclude='qiime2/**/_tmp_*' \
#   --exclude='qiime2/**/_export*' \

echo "[INFO] Done at: $(date)"
# echo "[INFO] Archive: qiime2_results.gz"
```
## 7) Outputs (high level)

`qiime2/demux.qzv` — read quality/demux report

`qiime2/table.qza` — ASV table (feature counts)

`qiime2/rep-seqs.qza` — ASV representative sequences

`qiime2/taxonomy/taxonomy_primary.qza` & `.qzv` — SILVA taxonomy + barplot

`qiime2/taxonomy/taxonomy_primary.tsv` — exported taxonomy mapping

`qiime2/feature-table.tsv` — exported ASV table (TSV)

**Summaries (SILVA):**

`qiime2/summary/asv_{phylum|family|genus|species}.tsv` (wide)

`qiime2/summary/asv_{rank}.long.tsv` (long)

`qiime2/summary/asv_{rank}.relabund.tsv` (relative abundance by sample)

**OTU97 (from ASV clustered at 97%)**:

`qiime2/otu97/otu97-table.qza`, `qiime2/otu97/otu97-rep-seqs.qza`

`qiime2/summary/otu97_{rank}.tsv` (+ `.long.tsv`, `.relabund.tsv`)

**GG2 counterparts** (if enabled): `asv_gg2_*`, `otu97_gg2_*`
**Concordance**: `qiime2/concordance/top50_concordance.csv`
**Phylogeny (optional)**: `qiime2/trees/gg2_nonv4_sepp/*`, including `insertion-tree.nwk`

## 8) Repro & troubleshooting

Rerun when rules/scripts change but inputs are fresh:
```
conda run -n snakemake snakemake --configfile workflow/config.yaml \
  --cores 8 --use-conda --printshellcmds \
  --rerun-triggers code  <targets...>

```
Clean specific outputs: delete files and re-invoke; or do a global clean with care:
```
conda run -n snakemake snakemake --configfile workflow/config.yaml \
  --delete-all-output --cores 1
```

Validate QIIME artifacts (`qiime tools validate file.qza`).

Ensure `metadata/sample-metadata.tsv` is **tab-separated** with first column `sample-id`.

**Primer consistency**: if primers change, **retrain** SILVA & GG2 classifiers.

## 9) Notes

**Primary taxonomy**: SILVA V3–V4 classifier trained from SSURef NR99 138.2 via RESCRIPt (extract → dereplicate uniq → train).

**Cross-check**: independent GG2 V3–V4 classifier built the same way.

**Rationale for** `symbiont_taxa.txt`: a curated endosymbiont list enabling (i) optional focused classifier and (ii) fast downstream subsetting of symbiont lineages for ecological analyses.

