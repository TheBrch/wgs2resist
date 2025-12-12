import pandas as pd
import os
from natsort import natsort_keygen
import argparse

parser = argparse.ArgumentParser(
    description="Generates an output file similar to snippy-core tab, just with non-synonymous mutations."
)
parser.add_argument("--samples", "-s", nargs="+", help="List of snippy tab files")
parser.add_argument("--outdir", "-o", help="Output directory")
args = parser.parse_args()

tables = []
dn_ds_set = []
annotations = pd.DataFrame()

for f in args.samples:
    sample = os.path.basename(os.path.dirname(f))
    tab = pd.read_csv(f)
    has_effect = tab[~tab.EFFECT.isna()]
    synonymous = has_effect.assign(
        SAMPLE=sample, SYN=has_effect.EFFECT.str.contains("synonymous_variant")
    )[["SAMPLE", "SYN", "LOCUS_TAG"]]
    dn_ds_set.append(synonymous)

    annot = has_effect[["LOCUS_TAG", "FTYPE", "GENE", "PRODUCT"]].drop_duplicates(
        "LOCUS_TAG"
    )
    if not annotations.empty:
        annot = annot[~annot["LOCUS_TAG"].isin(annotations["LOCUS_TAG"])]
    annotations = pd.concat([annotations, annot], ignore_index=True)

    no_synonymous = has_effect[~synonymous.SYN]
    no_complex = no_synonymous[~no_synonymous.TYPE.isin(["complex", "mnp", "ins"])]
    snd = no_complex[no_complex.REF.str.len() <= 2]

    mask = snd["REF"].str.len() > 1
    snd.loc[mask, "REF"] = snd.loc[mask, "REF"].str[-1]
    snd.loc[mask, "POS"] = snd.loc[mask, "POS"] + 1
    snd.loc[mask, "ALT"] = "-"

    simplified = snd.assign(
        CHR=snd.CHROM, SAMPLE=sample, SYN=snd.EFFECT.str.contains("synonymous_variant")
    )[["SAMPLE", "CHR", "POS", "REF", "ALT", "LOCUS_TAG"]]
    tables.append(simplified)


dn_ds = pd.concat(dn_ds_set, ignore_index=True)

dn_ds_summary = (
    dn_ds.groupby("LOCUS_TAG")["SYN"]
    .agg(
        n_synonymous=lambda x: x.sum(),
        n_nonsynonymous=lambda x: (~x).sum(),
        total_mutations="count",
    )
    .assign(dN_dS=lambda d: d["n_nonsynonymous"] / d["n_synonymous"].replace(0, pd.NA))
    .reset_index()
)
dn_ds_summary = dn_ds_summary.merge(annotations, on="LOCUS_TAG", how="left")

dn_ds_summary.to_csv(os.path.join(args.outdir, "dn_ds.tsv"), sep="\t", index=False)

collection = pd.concat(tables, ignore_index=True)

group_counts = collection.groupby(["CHR", "POS"]).size().reset_index(name="count")
valid_groups = group_counts[group_counts["count"] >= 0.05 * len(tables)]

filtered_collection = collection.merge(
    valid_groups[["CHR", "POS"]], on=["CHR", "POS"], how="inner"
)

wide_output = filtered_collection.pivot_table(
    index=["CHR", "POS", "REF"], columns="SAMPLE", values="ALT", aggfunc="first"
).reset_index()


sample_columns = [
    col for col in wide_output.columns if col not in ["CHR", "POS", "REF"]
]

for col in sample_columns:
    wide_output[col] = wide_output[col].fillna(wide_output["REF"])

wide_output = wide_output.sort_values(
    by=["CHR", "POS"], key=natsort_keygen()
).reset_index(drop=True)


wide_output.to_csv(
    os.path.join(args.outdir, "non-synonymous-core.tsv"), sep="\t", index=False
)
