#commiting to git create new project under version control link to already excisting file on desktop link to git repo url,
#developer settings generate new token selceted scope
#used prepossed data ie primer sequence removed
library(dada2)
list.files("~/Desktop/DADA2 metagenomics - basic SS") #check exact file name then copy it
path <- "~/Desktop/DADA2 metagenomics - basic SS/MiSeq_SOP-2"
list.files(path)

#matched list for forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


#retaining sample names. cutting the nmes at each underscore(since its sepatred by them)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #c(1,3) can be sued to select more than one section corresponding to the placement
sample.names <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

#qualityPLOTS
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#ensuring proper overlap- your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.

#new file with sample names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
#truncLen-toss reads with bp more than x with F and R respectively
#maxN-removes any N bases dad2 cant handle it
#rm.phix to remove the control genome virus of illumina
#maxEE-remove error rate above 2 for each read
#truncQ- read from begining and when error rate of 2 or worse is met toos anything after that, iff too short reead is tossed.

head(out)

#learn error rates
#eroor models for F and R
#algo using a rough guess ie common as real and variants as non real checks the sequences again and determines which is real, repeats
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#applying errorrates to data 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err =errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

#Merging reads
#OVERlap of at lest 12 bases or justConcatentate=TRUE
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#Amplicon sequencing variant
seqtab <- makeSequenceTable(mergers)
dim(seqtab)#sequencingtable
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)#sequencing table no chimera

sum(seqtab.nochim)/sum(seqtab)

#checking prpgress so far 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


#Assining taxonmy
taxa <- assignTaxonomy(seqtab.nochim,"~/Desktop/DADA2 metagenomics - basic SS/silva_nr99_v138.2_toGenus_trainset.fa", multithread = TRUE )
taxa.print <- taxa
rownames(taxa.print) <- NULL #removing seq rows for visuals
head(taxa.print)


###.    P H Y L O S E Q VISULIZATION   ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library(phyloseq)
library(Biostrings)
library(ggplot2)

#Df to use
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

#phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

#converting to short names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
