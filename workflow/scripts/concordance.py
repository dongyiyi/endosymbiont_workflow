import pandas as pd

table_tsv = snakemake.input["table_tsv"]
silva_tsv = snakemake.input["silva_tsv"]
gg2_tsv   = snakemake.input["gg2_tsv"]
out_csv   = snakemake.output["csv"]

counts = pd.read_csv(table_tsv, sep='\t', skiprows=1)
counts.rename(columns={counts.columns[0]:"Feature ID"}, inplace=True)
count_sum = counts.set_index("Feature ID").sum(axis=1).sort_values(ascending=False)

def load_tax(path):
    df = pd.read_csv(path, sep='\t')
    keep = [c for c in ["Feature ID","Taxon","Confidence"] if c in df.columns]
    return df[keep]

silva = load_tax(silva_tsv).set_index("Feature ID")
gg2   = load_tax(gg2_tsv).set_index("Feature ID")

TOP_N =  int(snakemake.params.get("top_n", 50)) if hasattr(snakemake, "params") else 50
top_ids = count_sum.index[:TOP_N]

out = pd.DataFrame({
    "Feature ID": top_ids,
    "Total Count": count_sum.loc[top_ids].values,
    "SILVA Taxon": silva.reindex(top_ids)["Taxon"].values,
    "SILVA Confidence": silva.reindex(top_ids)["Confidence"].values,
    "GG2 Taxon": gg2.reindex(top_ids)["Taxon"].values,
    "GG2 Confidence": gg2.reindex(top_ids)["Confidence"].values,
})

def rank_path(taxon):
    if pd.isna(taxon): return []
    return [t.strip() for t in str(taxon).split(";")]

def agreement_level(a, b):
    pa, pb = rank_path(a), rank_path(b)
    if not pa or not pb: return "none"
    levels = ["kingdom","phylum","class","order","family","genus","species"]
    same = 0
    for i in range(min(len(pa), len(pb))):
        if pa[i] == pb[i] and pa[i] != "Unassigned":
            same = i+1
        else:
            break
    return "none" if same == 0 else levels[min(same-1, len(levels)-1)]

out["Agreement"] = [agreement_level(a, b) for a,b in zip(out["SILVA Taxon"], out["GG2 Taxon"])]
out.to_csv(out_csv, index=False)
