# endosymbiont workflow
Endosymbiont 16S V3–V4 Workflow (Snakemake + QIIME 2 [https://qiime2.org/])
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
set -euo pipefail

CONFIG="${CONFIG:-workflow/config.yaml}"
CORES="${CORES:-8}"
WITH_GG2="${WITH_GG2:-no}"     # yes/no
WITH_TREE="${WITH_TREE:-no}"   # yes/no
EXTRA_ARGS="${EXTRA_ARGS:-}"   # e.g. "--rerun-incomplete --keep-going"

CORE_TARGETS=(
  qiime2/demux.qzv
  qiime2/table.qza
  qiime2/rep-seqs.qza
  qiime2/taxonomy/taxa-barplot_primary.qzv
  qiime2/taxonomy/taxonomy_primary.tsv
  qiime2/feature-table.tsv

  qiime2/summary/asv_genus.tsv
  qiime2/summary/asv_genus.long.tsv
  qiime2/summary/asv_genus.relabund.tsv
  qiime2/summary/asv_species.tsv
  qiime2/summary/asv_species.long.tsv
  qiime2/summary/asv_species.relabund.tsv

  qiime2/otu97/otu97-table.qza
  qiime2/otu97/otu97-rep-seqs.qza
  qiime2/summary/otu97_genus.tsv
  qiime2/summary/otu97_genus.long.tsv
  qiime2/summary/otu97_genus.relabund.tsv
  qiime2/summary/otu97_species.tsv
  qiime2/summary/otu97_species.long.tsv
  qiime2/summary/otu97_species.relabund.tsv

  qiime2/concordance/top50_concordance.csv
)

GG2_TARGETS=(
  qiime2/taxonomy/taxonomy_gg2.qza
  qiime2/taxonomy/taxonomy_gg2.tsv

  qiime2/summary/asv_gg2_genus.tsv
  qiime2/summary/asv_gg2_genus.long.tsv
  qiime2/summary/asv_gg2_genus.relabund.tsv
  qiime2/summary/asv_gg2_species.tsv
  qiime2/summary/asv_gg2_species.long.tsv
  qiime2/summary/asv_gg2_species.relabund.tsv

  qiime2/otu97/taxonomy_otu97_gg2.qza
  qiime2/summary/otu97_gg2_genus.tsv
  qiime2/summary/otu97_gg2_genus.long.tsv
  qiime2/summary/otu97_gg2_genus.relabund.tsv
  qiime2/summary/otu97_gg2_species.tsv
  qiime2/summary/otu97_gg2_species.long.tsv
  qiime2/summary/otu97_gg2_species.relabund.tsv
)

TREE_TARGETS=(
  qiime2/trees/gg2_nonv4_sepp/insertion-tree.qza
  qiime2/trees/gg2_nonv4_sepp/placements.qza
  qiime2/trees/gg2_nonv4_sepp/insertion-tree.nwk
)

TARGETS=("${CORE_TARGETS[@]}")
[[ "${WITH_GG2}" == "yes" ]]  && TARGETS+=("${GG2_TARGETS[@]}")
[[ "${WITH_TREE}" == "yes" ]] && TARGETS+=("${TREE_TARGETS[@]}")

echo "[run_all] config=${CONFIG} cores=${CORES} gg2=${WITH_GG2} tree=${WITH_TREE}"
conda run -n snakemake snakemake \
  --configfile "${CONFIG}" \
  --cores "${CORES}" \
  --use-conda --printshellcmds ${EXTRA_ARGS} \
  "${TARGETS[@]}"
```

Usage:
```
# SILVA only
./run_all.sh

# + GG2 tables
WITH_GG2=yes ./run_all.sh

# + GG2 + phylogeny
WITH_GG2=yes WITH_TREE=yes ./run_all.sh

# Resume robustly after interruption
EXTRA_ARGS="--rerun-incomplete --keep-going" ./run_all.sh
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

**Consensus handling** (optional, if enabled in your Snakefile): we provide per-DB tables and a concordance summary; a conservative consensus table can be produced that marks disagreements as `Unresolved@rank` (or applies genus-level fallback policies if configured).
