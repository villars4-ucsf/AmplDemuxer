.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0-CBI")
library(dada2)
library(tidyverse)

temp1 <- paste0("PARAV3_Val1n2_", c("1A","1B","2"), "/results/demux/")
temp2 <- paste0(list.files(path = "PARAV3_Val1n2_1B/results/demux",pattern = "S."),"/trimmed_noprimers")
temp3 <- c(paste0(temp1[1],temp2),paste0(temp1[2],temp2),paste0(temp1[3],temp2))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path=temp3, pattern="_R1_trimmed_noprimers.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path=temp3, pattern="_R2_trimmed_noprimers.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

filtFs <- paste0(sapply(strsplit(fnFs,"/P"),"[",1),"/filtered/",sample.names,"_F_filt.fastq.gz")
filtRs <- paste0(sapply(strsplit(fnRs,"/P"),"[",1),"/filtered/",sample.names,"_R_filt.fastq.gz")

names(filtFs) <- sample.names
names(filtRs) <- sample.names


# fnFs is path to reads
###########################################
# 1. filter and trim all files
out <- filterAndTrim(fnFs, #input
                     filtFs, #filtering output
                     fnRs,
                     filtRs,
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=c(5,5),
                     rm.phix=TRUE,
                     compress=TRUE,
                     multithread=TRUE,
                     trimRight = c(0,0),trimLeft = 1, minLen=75,matchIDs=TRUE) 
head(out)

# 2. learn errors using all files
#    goes 
errF <- learnErrors(filtFs[out[,2]>0], #remove files that have zero reads after filtering
                    multithread=TRUE,
                    MAX_CONSIST=10,
                    randomize=TRUE)
errR <- learnErrors(filtRs[out[,2]>0],
                    multithread=TRUE,
                    MAX_CONSIST=10,
                    randomize=TRUE)

derepFs <- derepFastq(filtFs[out[,2]>0],
                      verbose = TRUE)
derepRs <- derepFastq(filtRs[out[,2]>0],
                      verbose = TRUE)

dadaFs <- dada(derepFs, err=errF, selfConsist=TRUE, multithread=TRUE, verbose=TRUE, OMEGA_A=1e-120)
dadaRs <- dada(derepRs, err=errR, selfConsist=TRUE, multithread=TRUE, verbose=TRUE, OMEGA_A=1e-120)



mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=FALSE, trimOverhang = TRUE,minOverlap=1)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
temp1=data.frame(a = sapply(dadaFs,getN),samples=names(sapply(dadaFs, getN)))
temp2=data.frame(a = sapply(dadaRs,getN),samples=names(sapply(dadaRs, getN)))
temp3=data.frame(a = sapply(mergers,getN),samples=names(sapply(mergers, getN)))
temp4=data.frame(a = rowSums(seqtab.nochim),samples=names(rowSums(seqtab.nochim)))

out2 = as.data.frame(out)
out2$samples = rownames(out)
out2$samples =sapply(strsplit(out2$samples,"_R1_trimmed_noprimers.fastq.gz"),"[",1)

track = left_join(out2 ,temp1,by="samples") %>% 
	left_join(temp2,by="samples") %>% 
	left_join(temp3,by="samples") %>% 
	left_join(temp4,by="samples") %>%
	select(-samples)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
save(out,dadaFs, dadaRs,mergers,seqtab,seqtab.nochim, file = "PARAV3_Val1n2.RData")