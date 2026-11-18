"""
prinseq.py - Quality filtering of reads with PRINSEQ-lite (Snakemake script)

Removes reads containing Ns (-ns_max_n 0) and reads below the mean quality
threshold from config["qc"]["mq"] (-min_qual_mean). Works for paired-end
(2 input files) and single-end (1 input file) data.

PRINSEQ only takes an output PREFIX and appends its own file endings,
so the ending is cut off the output path expected by Snakemake.
"""

from subprocess import run

# --- Paired-end or single-end? ---------------------------------------------
# 2 input files (forward R1 + reverse R2) -> paired-end
# 1 input file                            -> single-end
paired = len(snakemake.input) == 2

# --- Output prefix -----------------------------------------------------------
# Cut off the file ending PRINSEQ appends itself:
#   paired-end: last 8 characters ("_1.fastq")  e.g. S1_A_1.fastq -> S1_A
#   single-end: last 6 characters (".fastq")    e.g. S1_A_1.fastq -> S1_A_1
cut = 8 if paired else 6
output_edit = str(snakemake.output[0])[:-cut]   # prefix for reads that PASS
output_bad = output_edit + "_low_qual"          # prefix for reads that FAIL

# --- Build PRINSEQ command ---------------------------------------------------
cmd = ["prinseq-lite.pl", "-verbose",            # print progress to stderr
       "-fastq", snakemake.input[0]]             # forward / single-end reads
if paired:
    cmd += ["-fastq2", snakemake.input[1]]       # reverse reads

cmd += ["-ns_max_n", "0",                        # QUALITY FILTER 1: no N allowed
        "-min_qual_mean", str(snakemake.params), # QUALITY FILTER 2: min. mean quality
        "-out_good", output_edit,                # prefix: passing reads
        "-out_bad", output_bad,                  # prefix: filtered-out reads
        "-log", str(snakemake.log)]              # PRINSEQ log file (read counts)

# --- Run PRINSEQ -------------------------------------------------------------
# check=True: stop the rule with an error if PRINSEQ fails
run(cmd, check=True)
