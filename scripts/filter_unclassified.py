import pandas as pd
import re

# Get input and output file paths from Snakemake
input_taxonomy = snakemake.input["raw_taxonomy"]
input_summary = snakemake.input["raw_summary"]
output_taxonomy = snakemake.output["cleaned_taxonomy"]
output_summary = snakemake.output["cleaned_summary"]

def clean_taxonomy(tax_string):
    # Handle missing values safely (e.g., NaN)
    if pd.isna(tax_string):
        return tax_string

    # Split taxonomy string into tokens; trim whitespace; drop empty tokens
    taxa = [t.strip() for t in str(tax_string).split(";") if t.strip()]

    cleaned = []
    first_unclassified = None

    for t in taxa:
        # Match tokens like "Bacillus_unclassified(57)" or "Bacillus_unclassified"
        if re.fullmatch(r".*_unclassified(\(\d+\))?", t):
            first_unclassified = t
            break  # Stop here: everything below is redundant/uncertain
        cleaned.append(t)

    # Keep exactly one "_unclassified" token (the first one) as an explicit marker
    if first_unclassified is not None:
        cleaned.append(first_unclassified)

    # Re-join the cleaned taxonomy tokens
    return ";".join(cleaned) + ";"

# Clean taxonomy file (no header)
df_tax = pd.read_csv(input_taxonomy, sep="\t", header=None)
df_tax[1] = df_tax[1].apply(clean_taxonomy)
df_tax.to_csv(output_taxonomy, sep="\t", index=False, header=False)

# Clean summary file (has header)
df_summary = pd.read_csv(input_summary, sep="\t")

# Apply the same cleaning to the 'taxonomy' column
df_summary["taxonomy"] = df_summary["taxonomy"].apply(clean_taxonomy)

df_summary.to_csv(output_summary, sep="\t", index=False)
