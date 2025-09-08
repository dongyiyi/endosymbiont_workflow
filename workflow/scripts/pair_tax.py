#!/usr/bin/env python
import argparse, pandas as pd, re, numpy as np

RANKS = ["kingdom","phylum","class","order","family","genus","species"]

def read_tax_tsv(p):
    df = pd.read_csv(p, sep="\t", dtype=str)
    if "Feature ID" in df.columns and "Taxon" in df.columns:
        df = df.rename(columns={"Feature ID":"FeatureID"})
    elif "FeatureID" not in df.columns:
        df = df.rename(columns={df.columns[0]:"FeatureID"})
    if "Taxon" not in df.columns:
        raise SystemExit(f"[ERR] {p} missing Taxon column")
    if "Confidence" not in df.columns:
        df["Confidence"] = ""
    df["Taxon"] = df["Taxon"].fillna("")
    return df[["FeatureID","Taxon","Confidence"]]

def split_taxon(t):
    parts = [x.strip() for x in re.split(r";\s*", t) if x.strip()]
    ranks = [""]*7
    for p in parts:
        if "__" in p:
            a,b = p.split("__",1)
            key = a.lower()
            if key.startswith("k") or key in ("d","d_0","k_0"): idx=0
            elif key.startswith("p") or key=="d_1": idx=1
            elif key.startswith("c") or key=="d_2": idx=2
            elif key.startswith("o") or key=="d_3": idx=3
            elif key.startswith("f") or key=="d_4": idx=4
            elif key.startswith("g") or key=="d_5": idx=5
            elif key.startswith("s") or key=="d_6": idx=6
            else: continue
            ranks[idx] = b
    return pd.Series(ranks, index=RANKS)

def deepest_agree(row):
    d = None
    for i,r in enumerate(RANKS):
        a = row[f"silva_{r}"]; b = row[f"gg2_{r}"]
        if a and b and a==b: d = r
        else: break
    return d if d else ""

def first_diff(row):
    for r in RANKS:
        a = row[f"silva_{r}"]; b = row[f"gg2_{r}"]
        if (a or b) and (a != b):
            return r
    return ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--silva", required=True)
    ap.add_argument("--gg2", required=True)
    ap.add_argument("--out-pairs", required=True)
    ap.add_argument("--out-agree", required=True)
    args = ap.parse_args()

    s = read_tax_tsv(args.silva)
    g = read_tax_tsv(args.gg2)

    s = s.join(s["Taxon"].apply(split_taxon))
    g = g.join(g["Taxon"].apply(split_taxon))
    s = s.add_prefix("silva_").rename(columns={"silva_FeatureID":"FeatureID"})
    g = g.add_prefix("gg2_").rename(columns={"gg2_FeatureID":"FeatureID"})

    df = s.merge(g, on="FeatureID", how="outer").fillna("")
    for r in RANKS:
        df[f"agree_{r}"] = (df[f"silva_{r}"]!="" ) & (df[f"gg2_{r}"]!="") & (df[f"silva_{r}"]==df[f"gg2_{r}"])

    df["deepest_agree"] = df.apply(deepest_agree, axis=1)
    df["first_diff"]    = df.apply(first_diff, axis=1)

    # 输出并列清单
    keep_cols = (["FeatureID","silva_Taxon","silva_Confidence","gg2_Taxon","gg2_Confidence"] +
                 sum(([f"silva_{r}",f"gg2_{r}",f"agree_{r}"] for r in RANKS), [] ) +
                 ["deepest_agree","first_diff"])
    df[keep_cols].to_csv(args.out_pairs, sep="\t", index=False)

    # 各阶元一致性统计
    rows=[]
    for r in RANKS:
        a = df[f"silva_{r}"]; b = df[f"gg2_{r}"]
        both_called = (a!="") & (b!="")
        n_both = int(both_called.sum())
        n_agree = int((df[f"agree_{r}"] & both_called).sum())
        rows.append({
            "rank": r,
            "n_total": int(len(df)),
            "n_both_called": n_both,
            "n_agree": n_agree,
            "n_conflict": n_both - n_agree,
            "n_silva_only": int(((a!="") & (b=="")).sum()),
            "n_gg2_only": int(((a=="") & (b!="")).sum()),
            "agree_rate_both_called": (n_agree/n_both) if n_both else np.nan
        })
    pd.DataFrame(rows).to_csv(args.out_agree, sep="\t", index=False)

if __name__ == "__main__":
    main()
