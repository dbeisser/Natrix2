library(dada2)
library(ShortRead)

# Set up log files
log <- file(toString(snakemake@log), open="wt")
sink(log, append = TRUE)
sink(log, type="message", append = TRUE)

# Get parameters from Snakemake
paired.end <- snakemake@params[["paired_end"]]
min.overlap <- snakemake@params[["minoverlap"]]
split.samples <- snakemake@params[["splitsamples"]] == "split_sample"
sample.names <- sapply(strsplit(basename(snakemake@input[["fwd"]]), "_[12]_cut.fastq"), `[`, 1)
fnFs <- snakemake@input[["fwd"]]
names(fnFs) <- sample.names

print(paired.end)
print(split.samples)
print(sample.names)

set.seed(100)

# Learn to guess errors
errF <- learnErrors(fnFs, nbases=1e8, multithread=TRUE, randomize=TRUE, verbose=TRUE)
if (paired.end) {
  fnRs <- snakemake@input[["rev"]]
  names(fnRs) <- sample.names
  errR <- learnErrors(fnRs, nbases=1e8, multithread=TRUE, randomize=TRUE, verbose=TRUE)
}

if (split.samples) {
  samples <- unlist(unique(strsplit(sample.names, "_[AB]")))
  sample.names2 <- lapply(samples, function(x) grep(x, sample.names, value=TRUE))
  names(sample.names2) <- samples
} else {
  sample.names2 <- sample.names
  names(sample.names2) <- sample.names
}

# Access the current sample
sam <- sample.names2[[1]]
result <- c()

sample.name <- snakemake@wildcards[["sample"]]
unit.name <- snakemake@wildcards[["unit"]]
cat("Processing:", sample.name, unit.name, "\n")

# Dereplication and denoising
dadaF <- dada(fnFs[unlist(sam)], err=errF, multithread=TRUE, verbose=TRUE)
print("Forward After DADA2 dereplication and denoising:")
print(dadaF)

if (paired.end) {
  dadaR <- dada(fnRs[unlist(sam)], err=errR, multithread=TRUE, verbose=TRUE)
  print("Reverse After DADA2 dereplication and denoising:")
  print(dadaR)
  
  result <- mergePairs(dadaF, fnFs[unlist(sam)], dadaR, fnRs[unlist(sam)], verbose=TRUE, minOverlap=min.overlap)
  if (!is(result, "list")) {
    result <- list(result)
  }
} else {
  if (is(dadaF, "dada")) {
    dadaF <- list(dadaF)
  }
  seqtab <- lapply(dadaF, function(x) t(makeSequenceTable(x)))
  result <- lapply(seqtab, function(x) data.frame(abundance=x[,1], sequence=rownames(x)))
}

print("Sequences left after assembly:")
for (x in result) {
  print(nrow(x))
}

for (j in seq_along(result)) {
  if (nrow(result[[j]]) == 0) {
    print(paste0("No sequences left for sample ", names(sam)))
    cat(NULL, file=snakemake@output[[j]])
  } else {
    seq_names <- paste0(seq(nrow(result[[j]])), ";size=", result[[j]]$abundance, ";")
    uniquesToFasta(result[[j]], fout=snakemake@output[[j]], ids=seq_names)
  }
}

sink()
sink(type="message")
