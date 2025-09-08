#!/usr/bin/env python
import argparse, pandas as pd, re

RANKS = ["kingdom","phylum","class","order","family","genus","species"]
PREFIX = ["k__","p__","c__","o__","f__","g__","s__"]

def read_tax(p):
    df = pd.read_csv(p, sep="\t", dtype=str)
    if "Feature ID" in df.columns and "Taxon" in df.columns:
        df = df.rename(columns={"Feature ID":"FeatureID"})
    elif "FeatureID" not in df.columns:
        df = df.rename(columns={df.columns[0]:"FeatureID"})
    if "Taxon" not in df.columns: raise SystemExit(f"{p} missing Taxon")
    if "Confidence" not in df.columns: df["Confidence"] = ""
    df["Taxon"] = df["Taxon"].fillna("")
    return df[["FeatureID","Taxon","Confidence"]]

def split_taxon(t):
    parts = [x.strip() for x in re.split(r";\s*", t) if x.strip()]
    ranks = [""]*7
    for p in parts:
        if "__" in p:
            a,b = p.split("__",1)
            key=a.lower()
            idx = {"k":0,"d_0":0,"p":1,"d_1":1,"c":2,"d_2":2,"o":3,"d_3":3,"f":4,"d_4":4,"g":5,"d_5":5,"s":6,"d_6":6}.get(key[0],None)
            if idx is None: 
                if key in ("d","k_0"): idx=0
                else: continue
            ranks[idx]=b
    return ranks

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--silva", required=True)
    ap.add_argument("--gg2", required=True)
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    s = read_tax(args.silva)
    g = read_tax(args.gg2)

    s_ranks = s["Taxon"].apply(split_taxon).to_list()
    g_ranks = g["Taxon"].apply(split_taxon).to_list()

    s = s.assign(**{r:[sr[i] for sr in s_ranks] for i,r in enumerate(RANKS)})
    g = g.assign(**{r:[gr[i] for gr in g_ranks] for i,r in enumerate(RANKS)})

    df = s.merge(g, on="FeatureID", how="outer", suffixes=("_silva","_gg2")).fillna("")

    out_rows=[]
    for _,row in df.iterrows():
        # 寻找“最深一致阶元”
        deepest = -1
        for i,r in enumerate(RANKS):
            a=row[f"{r}_silva"]; b=row[f"{r}_gg2"]
            if a and b and a==b:
                deepest=i
            else:
                break
        if deepest>=0:
            lineage = [f"{PREFIX[i]}{row[f'{RANKS[i]}_silva']}" for i in range(deepest+1)]
            taxon = "; ".join(lineage)
            conf = max(row["Confidence_silva"], row["Confidence_gg2"]) if (row["Confidence_silva"] and row["Confidence_gg2"]) else (row["Confidence_silva"] or row["Confidence_gg2"])
        else:
            # 没有任何阶元达成一致：保守置为未定
            taxon = ""
            conf = row["Confidence_silva"] or row["Confidence_gg2"]
        out_rows.append({"FeatureID":row["FeatureID"], "Taxon":taxon, "Confidence":conf})

    pd.DataFrame(out_rows)[["FeatureID","Taxon","Confidence"]].to_csv(args.out, sep="\t", index=False)

if __name__=="__main__":
    main()
