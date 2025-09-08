#!/usr/bin/env bash
set -euo pipefail

# 1) mini raw fastq
mkdir -p tests/data/raw
cat > tests/data/Sample1_R1.fastq <<'EOF'
@S1/1
ACGTACGTACGT
+
FFFFFFFFFFFF
EOF
cat > tests/data/Sample1_R2.fastq <<'EOF'
@S1/2
ACGTACGTACGT
+
FFFFFFFFFFFF
EOF
gzip -f tests/data/Sample1_R1.fastq
gzip -f tests/data/Sample1_R2.fastq
mkdir -p raw
cp -f tests/data/Sample1_R1.fastq.gz raw/Sample1_R1.fastq.gz
cp -f tests/data/Sample1_R2.fastq.gz raw/Sample1_R2.fastq.gz

# 2) mini metadata（QIIME2：sample-id；group）
mkdir -p metadata
cat > metadata/sample-metadata.tsv <<'EOF'
sample-id\tgroup
Sample1\tnative
EOF

# 3) classifier / refs（only for dry-run）
mkdir -p refs/V3V4-silva-138.2 refs/GG2-V3V4 refs/GG2-raw
touch refs/V3V4-silva-138.2/silva-138.2-V3V4-uniq-classifier.qza
touch refs/GG2-V3V4/V3V4-uniq-classifier.qza
# 
touch refs/GG2-raw/2024.09.backbone.full-length.fna.qza
touch refs/GG2-raw/2024.09.backbone.tax.qza
touch refs/GG2-raw/2022.10.backbone.sepp-reference.qza

# 4) CI 
mkdir -p tests
cat > tests/config.ci.yaml <<'EOF'
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

# 
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
  # 
  classifier_qza: refs/GG2-V3V4/V3V4-uniq-classifier.qza
  outdir: refs/GG2-V3V4
  seqs_qza: refs/GG2-raw/2024.09.backbone.full-length.fna.qza
  tax_qza:  refs/GG2-raw/2024.09.backbone.tax.qza
  sepp_ref_qza: refs/GG2-raw/2022.10.backbone.sepp-reference.qza

summary:
  group_column: group
concordance:
  top_n: 10
EOF
