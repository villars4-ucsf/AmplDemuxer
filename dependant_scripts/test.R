

library(dada2)

# path to data
path <- "~/eppicenter/library/WDL/AmplDemuxer/dependant_scripts/scratch/MiSeq_SOP/"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs,
                     filtFs,
                     fnRs,
                     filtRs,
                     truncLen=c(240,160),
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=2,
                     rm.phix=TRUE,
                     compress=TRUE,
                     multithread=TRUE) # On Windows set multithread=FALSE
