#!/usr/bin/env python
import argparse, pandas as pd
from pathlib import Path

TAX_COLS = {"Taxon","Kingdom","Phylum","Class","Order","Family","Genus","Species"}

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wide-tsv", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--group-col", required=True, help="metadata 分组列名，如 group / region / native_status")
    ap.add_argument("--out-rel",  required=True)
    ap.add_argument("--out-long", required=True)
    return ap.parse_args()

def main():
    a = parse_args()
    out_rel  = Path(a.out_rel);  out_rel.parent.mkdir(parents=True, exist_ok=True)
    out_long = Path(a.out_long); out_long.parent.mkdir(parents=True, exist_ok=True)

    wide = pd.read_csv(a.wide_tsv, sep="\t")
    meta = pd.read_csv(a.metadata, sep="\t", dtype=str)

    # 识别样本列
    sample_cols = [c for c in wide.columns if c not in TAX_COLS]
    assert sample_cols, "未找到样本列，请检查宽表列名。"

    # 确认 metadata 有 sample ID 列（QIIME2 常用 'sample-id' / 'SampleID' / '#SampleID'）
    id_candidates = ["sample-id","SampleID","#SampleID","#Sample ID","sample_name","id","ID"]
    id_col = next((c for c in id_candidates if c in meta.columns), None)
    if id_col is None:
        raise SystemExit(f"metadata 缺少样本ID列（期望之一：{id_candidates}）。")

    if a.group_col not in meta.columns:
        raise SystemExit(f"metadata 缺少分组列 {a.group_col}。")

    # 计算相对丰度（按列每个样本除以该样本总和）
    rel = wide.copy()
    sums = rel[sample_cols].sum(axis=0).replace(0, 1)
    rel[sample_cols] = rel[sample_cols].div(sums, axis=1)
    rel.to_csv(out_rel, sep="\t", index=False)

    # 生成长表并附上分组
    long = rel.melt(id_vars=[c for c in rel.columns if c not in sample_cols],
                    value_vars=sample_cols,
                    var_name="SampleID", value_name="RelAbundance")
    meta_sub = meta[[id_col, a.group_col]].rename(columns={id_col:"SampleID"})
    long = long.merge(meta_sub, on="SampleID", how="left")
    long.to_csv(out_long, sep="\t", index=False)
    print(f"[wide_post] wrote {out_rel} and {out_long}")

if __name__ == "__main__":
    main()
