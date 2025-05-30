import pandas as pd
import re

# Get input and output file paths from Snakemake
input_taxonomy = snakemake.input["raw_taxonomy"]
input_summary = snakemake.input["raw_summary"]
output_taxonomy = snakemake.output["cleaned_taxonomy"]
output_summary = snakemake.output["cleaned_summary"]

# Clean taxonomy file (no header)
df_tax = pd.read_csv(input_taxonomy, sep="\t", header=None)

def clean_taxonomy(tax_string):
    taxa = tax_string.split(";")
    cleaned = [t for t in taxa if not re.match(r".*_unclassified(\(\d+\))?$", t.strip())]
    return ";".join(cleaned)

df_tax[1] = df_tax[1].apply(clean_taxonomy)
df_tax.to_csv(output_taxonomy, sep="\t", index=False, header=False)

# Clean summary file (has header)
df_summary = pd.read_csv(input_summary, sep="\t")

# Apply the same cleaning to the 'taxonomy' column
df_summary["taxonomy"] = df_summary["taxonomy"].apply(clean_taxonomy)

df_summary.to_csv(output_summary, sep="\t", index=False)
