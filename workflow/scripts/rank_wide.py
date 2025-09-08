#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import re
from pathlib import Path

RANKS = ["kingdom","phylum","class","order","family","genus","species"]

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-tsv", required=True)
    ap.add_argument("--tax-tsv",   required=True)
    ap.add_argument("--rank",      required=True, choices=["phylum","family","genus","species"])
    ap.add_argument("--output",    required=True)
    return ap.parse_args()

def split_taxon(t):
    if pd.isna(t):
        return pd.Series([""]*7, index=RANKS)
    parts = [x.strip() for x in re.split(r";\s*", str(t)) if x.strip()]
    ranks = [""]*7
    for p in parts:
        if "__" in p:
            a,b = p.split("__",1)
            key = a.lower()
            # 宽松匹配多种前缀（SILVA/GG2 兼容）
            if   key in ["k","d","d_0","k_0"]: idx=0
            elif key in ["p","d_1","p_1"]:    idx=1
            elif key in ["c","d_2","c_2"]:    idx=2
            elif key in ["o","d_3","o_3"]:    idx=3
            elif key in ["f","d_4","f_4"]:    idx=4
            elif key in ["g","d_5","g_5"]:    idx=5
            elif key in ["s","d_6","s_6"]:    idx=6
            else: continue
            ranks[idx] = b
    return pd.Series(ranks, index=RANKS)

def main():
    args = parse_args()
    out_p = Path(args.output)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    # 表格（biom->tsv 导出）
    tab = pd.read_csv(args.table_tsv, sep="\t", dtype=str, header=0, skiprows=1)
    # 统一第一列列名
    first_col = tab.columns[0].strip()
    if first_col.lstrip().startswith("#"):
        first_col = first_col.lstrip("#").strip()
    if first_col.lower().replace(" ", "") in {"otuid", "featureid", "otu_id", "feature_id"}:
        tab = tab.rename(columns={tab.columns[0]: "FeatureID"})
    else:
        tab = tab.rename(columns={tab.columns[0]: "FeatureID"})

    # 样本列并转为数值
    sample_cols = [c for c in tab.columns if c != "FeatureID"]
    tab[sample_cols] = tab[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # taxonomy TSV
    tax = pd.read_csv(args.tax_tsv, sep="\t", dtype=str)
    if "FeatureID" not in tax.columns:
        tax = tax.rename(columns={tax.columns[0]:"FeatureID"})
    if "Taxon" not in tax.columns:
        # 保险：某些导出格式首列即为 Taxon
        maybe = [c for c in tax.columns if c.lower().startswith("tax")]
        if maybe:
            tax = tax.rename(columns={maybe[0]:"Taxon"})
    tax["Taxon"] = tax["Taxon"].fillna("")

    tx = tax.join(tax["Taxon"].apply(split_taxon))

    # 合并
    df = tab.merge(tx[["FeatureID","Taxon"]+RANKS], on="FeatureID", how="left").fillna("")
    level_name = args.rank  # e.g. genus
    key = df[level_name].replace({"":"Unassigned","__":"Unassigned"})
    df = df.assign(Key=key)

    # 按样本求和
    reads = df[sample_cols].sum(axis=1)
    df["_w"] = reads
    agg_counts = df.groupby("Key", dropna=False)[sample_cols].sum()

    # 组装谱系信息（每个 Key 选权重最大的一条作为代表）
    meta_rows = []
    for k, sub in df.groupby("Key", dropna=False):
        if k == "Unassigned":
            lineage = {r:"" for r in RANKS}
            lineage[level_name] = "Unassigned"
        else:
            lineage = {}
            cand = sub[[*_ for _ in RANKS] + ["_w"]]
            for r in RANKS:
                vals = cand.groupby(r)["_w"].sum().sort_values(ascending=False)
                lineage[r] = vals.index[0] if len(vals) else ""
            lineage[level_name] = k
        row = {"Taxon": k}
        row.update({r.capitalize(): lineage[r] for r in RANKS})
        meta_rows.append(row)
    meta = pd.DataFrame(meta_rows).set_index("Taxon")

    out = meta.join(agg_counts, how="left").reset_index()
    cols = ["Taxon"] + [r.capitalize() for r in RANKS] + sample_cols
    out = out[cols]
    out.to_csv(out_p, sep="\t", index=False)
    print(f"[rank_wide] wrote {out_p} ({len(out)} rows, {len(sample_cols)} samples)")

if __name__ == "__main__":
    main()
