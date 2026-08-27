import re
import dinopy
import logging
import subprocess
import numpy as np
import pandas as pd

# Process sequencing reads depending on whether the input is paired-end or single-end.
# Paired-end reads are assembled with PANDAseq.
# Single-end reads are processed directly in Python: primer/barcode/polyN regions
# are removed and the remaining sequence length is checked against the configured limits.

# Load the primer table and use the Probe column as the lookup key for each sample.
primer_table = pd.read_csv(snakemake.input.primer_t, index_col="Probe",
        na_filter=False).to_dict("index")

# Paired-end processing: delegate read assembly and optional primer removal to PANDAseq.
if snakemake.params.paired_end:
    # If primers were already removed, run PANDAseq without primer sequences.
    if snakemake.params.prim_rm:
        subprocess.call(["pandaseq",
            "-f",snakemake.input[0], "-r", snakemake.input[1], "-B", "-a", "-F",
            "-g", str(snakemake.log),
            "-w", str(snakemake.output), "-N",
            "-T", str(snakemake.threads),
            "-t", str(snakemake.params.threshold),
            "-o", str(snakemake.params.minoverlap),
            "-l", str(snakemake.params.minlen),
            "-L", str(snakemake.params.maxlen),
            "-C" "min_phred:" + str(snakemake.params.minqual)])
    else:
        # If primers are still present, obtain forward and reverse primers for this sample.
        r1_primer = primer_table[snakemake.wildcards.sample + "_"
            + snakemake.wildcards.unit]["specific_forward_primer"]
        r2_primer = primer_table[snakemake.wildcards.sample + "_"
            + snakemake.wildcards.unit]["specific_reverse_primer"]

        subprocess.call(["pandaseq",
            "-f", snakemake.input[0], "-r", snakemake.input[1], "-B", "-a", "-F",
            "-g", str(snakemake.log),
            "-w", str(snakemake.output),"-N",
            "-p", r1_primer, "-q", r2_primer,
            "-T", str(snakemake.threads),
            "-t", str(snakemake.params.threshold),
            "-o", str(snakemake.params.minoverlap),
            "-l", str(snakemake.params.minlen),
            "-L", str(snakemake.params.maxlen),
            "-C" "min_phred:" + str(snakemake.params.minqual)])
else:
    # Single-end processing: trim primer/barcode/polyN regions and filter reads by length.
    logging.basicConfig(filename=str(snakemake.log),
            level=logging.DEBUG)
    iupac_dict_regex = {"M":"[AC]", "R":"[AG]", "W":"[AT]", "S":"[CG]",
            "Y":"[CT]","K":"[GT]", "V":"[ACG]", "H":"[ACT]",
            "D":"[AGT]", "B":"[CGT]", "X":"[ACGT]", "N":"[ACGT]"}

    # Process all reads in one single-end FASTQ file.
    # Reads with a matching primer/barcode and a valid post-trimming length are written
    # to the assembled FASTQ. All other reads are written unchanged to filtered_out.
    # Sequence and quality values must always be trimmed by the same offset to keep
    # the resulting FASTQ records valid.
    def primer_len_filter(path, sample):
        sequence = dinopy.FastqReader(path)
        assembled = dinopy.FastqWriter(path.rsplit("_", 1)[0]
                + "_assembled.fastq")
        filt_out = dinopy.FastqWriter(path.rsplit("_", 1)[0]
                + "_filtered_out.fastq")
        assembled_counter = 0
        filt_out_counter = 0
        prim_counter = 0
        assembled.open()
        filt_out.open()
        # Read sequence and quality information together so both remain synchronized.
        for read in sequence.reads(quality_values=True):
            name = read.name.decode()
            # trim_end is the number of leading bases removed from the sequence.
            match, seq, trim_end = check_for_match(read.sequence.decode(), sample)
            if match and snakemake.params.maxlen >= len(seq) >= snakemake.params.minlen:
                # Apply the identical trim offset to quality scores; otherwise the FASTQ
                # sequence and quality lines would have different lengths.
                assembled.write(seq.encode(), name.split(" ")[0].encode(), read.quality[trim_end:])
                assembled_counter += 1
            else:
                if not match:
                    prim_counter +=1
                filt_out.write(read.sequence, name.split(" ")[0].encode(), read.quality)
                filt_out_counter += 1
        logging.info("{}: {} sequences were kept, \
                {} sequences were filtered out. \
                {} sequences without primer or barcode matches were detected. \
                If too many, check primer table or set all_primers to TRUE".format(sample,
                    assembled_counter, filt_out_counter, prim_counter))
        assembled.close()
        filt_out.close()

    # Convert IUPAC ambiguity codes in barcode sequences to regular-expression patterns.
    def iupac_replace(sequence, iupac_dict):
        for i, j in iupac_dict_regex.items():
            sequence = sequence.replace(i, j)
        return sequence

    # Search for the expected primer/barcode combination and remove the complete
    # leading polyN + primer + barcode region when a match is found.
    # Test the expected position as well as shifts of +/-1 and +/-2 bases to account
    # for small inaccuracies in the polyN length specified in the primer table.
    # Return: (match_found, trimmed_sequence, trim_offset).
    def check_for_match(sequence, sample):
        if snakemake.params.prim_rm:
            return (True, sequence, 0)
        else:
            poly_prim_bar = [primer_table[sample][key] for key
                    in primer_table[sample].keys() if key in ["poly_N",
                        "specific_forward_primer", "Barcode_forward"]]
            prim_bar = re.compile(poly_prim_bar[1]
                    + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
            # Do not reject a read until every allowed positional offset was tested.
            for i in [0, 1, -1, 2, -2]:
                start = np.clip(len(primer_table[sample]["poly_N"]) + i,
                        a_min=0, a_max=None)
                end = np.clip(len("".join(poly_prim_bar)) + i, a_min=0,
                        a_max=None)
                if prim_bar.match(sequence[start : end]):
                    # Slice at the exact end position instead of using str.replace(),
                    # which could remove an identical sequence occurring later in the read.
                    return (True, sequence[end:], end)
            # No primer/barcode match was found at any allowed offset.
            return (False, sequence, 0)

    primer_len_filter(snakemake.input[0],
            snakemake.input[0].split("/")[-1].rsplit("_", 1)[0])
