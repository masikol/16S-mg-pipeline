#!/usr/bin/env Rscript --vanilla
# Script performs filtering, denoising, merging and taxonomic classification of PE reads.
# It's main output is a phyloseq object.
# Usage:
#   Rscript --vanilla filter-merge-classify.R INPUT_DIR METADATA_FILE

# == Get and check arguments ==
argsNum <- 2 # there must be 2 arguments
if (length(commandArgs(T)) != argsNum) {
  message(sprintf("Invalid number of arguments: %d (must be %d).",
                  length(commandArgs(T)), argsNum))
  message("Your aruments: ", appendLF = F)
  for (arg in commandArgs(T)) {
    message(sprintf("`%s` ", arg), appendLF = F)
  }
  message("\n")
  quit(save = "no", status = 1)
}
rm(argsNum)

# == Store command line argument to variables ==
# Get and check path to input dir (it's directory with fastq files after preprocess16S)
input.dir <- normalizePath(
  trimws(commandArgs(T)[1], which="right", whitespace = "/")
)
if (!dir.exists(input.dir)) {
  message(sprintf("Input directory does not exist: `%s`!", input.dir))
  quit(save = "no", status = 1)
}

# Set workdir to parent  dir of input directory
wrkdir <- dirname(input.dir)

# Get and check path to metadata file
metadata.fpath <- normalizePath(commandArgs(T)[2])
if (!file.exists(metadata.fpath)) {
  message(sprintf("Metadata file does not exist: `%s`!", metadata.fpath))
  quit(save = "no", status = 1)
}

# Print some stuff
cat(sprintf("Input directory: `%s`.\n", input.dir))
cat(sprintf("Metadata file: `%s`.\n", metadata.fpath))

# == Read metadata file ==
metadata.df <- read.csv(metadata.fpath, header = TRUE)
rownames(metadata.df) <- metadata.df$SampleID

# == Go to work dir ==
setwd(wrkdir)

# == Create directories for plots and output files ==
plots.dir <- file.path(wrkdir, "plots")
dir.create(plots.dir)
outfiles.dir <- file.path(wrkdir, "output_files")
dir.create(outfiles.dir)

# == Check if reference file for classification exists ==
ref.tax.file <- Sys.getenv("REF_CLASSIF_PATH")
if (! file.exists(ref.tax.file)) {
  message(sprintf("Reference file for taxonomic classification does not exist: `%s`!", ref.tax.file))
  quit(save = "no", status = 1)
} else {
  cat(sprintf("\nFound reference file for taxonomic classification: `%s`\n", ref.tax.file))
}

# == Import necessary packages ==
cat("\n")
cat("Importing phyloseq...\n")
library("phyloseq")
cat("Importing ggplot2...\n")
library("ggplot2")
cat("Importing dada2...\n")
library("dada2")
cat("Importing msa...\n")
library("msa")
cat("Importing phangorn...\n")
library("phangorn")

# == Get and check read files ==
fnFs <- sort(list.files(input.dir, pattern="_R1.*.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(input.dir, pattern="_R2.*.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

if (length(fnFs) != length(fnRs)) {
  message(sprintf("Nubmer of R1 files (%d) is not equal to number of R2 files (%d).",
                  length(fnFs), length(fnRs)))
  quit(save = "no", status = 1)
} else {
  cat("\n")
  cat(sprintf("%d files of forward reads detected.\n", length(fnFs)))
  cat(sprintf("%d files of reverse reads detected.\n", length(fnRs)))
}

# Check if all samples have their fastq files
if (!all(sample.names %in% metadata.df$SampleID)) {
  message('\n----------')
  message("Error: some samples have their fastq files but are not mentioned in metadata file!\n")
  message("Here are missing samples:\n")
  message(sample.names[ ! sample.names %in% metadata.df$SampleID])
  message('----------')
  quit(save = "no", status = 1)
}

# Check if all fasta files are mentioned in metadata file
if (!all(metadata.df$SampleID %in% sample.names)) {
  message('\n----------')
  message("Error: some samples mentioned in metadata file have no input fastq files!\n")
  message("Here are missing samples:\n")
  message(metadata.df$SampleID[ ! metadata.df$SampleID %in% sample.names])
  message('----------')
  quit(save = "no", status = 1)
}

# == Plot initial quality profiles ==
qual_plots_fpath <- file.path(plots.dir, "initial_quality_profiles.pdf")
cat("\n")
cat(sprintf("Saving initial quality profiles to `%s`\n", qual_plots_fpath))
cat(sprintf("0/%d done", length(fnFs)))

# Make plots and save them
pdf(qual_plots_fpath)
for (i in 1:length(fnFs)) {
  print(plotQualityProfile(fnFs[i]))
  print(plotQualityProfile(fnRs[i]))
  cat(sprintf("\r%d/%d done", i, length(fnFs)))
}
cat("\n")
dev.off()
cat("Inittial quality profiles are created.\n")

# == Filter and trim ==
cat("\n")
cat("Reads filtering started.\n")
filtFs <- file.path(wrkdir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(wrkdir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN = 0, maxEE = c(2.0, 2.0), rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE, verbose = TRUE)
cat("Reads filtering is completed.\n")

# == Plot quality profiles after QC ==
qual_plots_fpath <- file.path(plots.dir, "quality_profiles_after_QC.pdf")
cat("\n")
cat(sprintf("Saving quality profiles after QC to `%s`\n", qual_plots_fpath))
cat(sprintf("0/%d done", length(fnFs)))

# Make plots ans save them
pdf(qual_plots_fpath)
for (i in 1:length(filtFs)) {
  print(plotQualityProfile(filtFs[i]))
  print(plotQualityProfile(filtRs[i]))
  cat(sprintf("\r%d/%d done", i, length(filtFs)))
}
cat("\n")
dev.off()
cat("Quality profiles after QC are created.\n")

# == Dereplicate ==
cat("\n")
cat("Dereplication started.\n")
derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.names
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepRs) <- sample.names
cat("Dereplication is completed.\n")

rm(filtFs, filtRs)

# == Learn the Error Rates ==
cat("\n")
cat("Learning error models started (it will take a while).\n")
errF <- learnErrors(derepFs, multithread = T, verbose = T)
errR <- learnErrors(derepRs, multithread = T, verbose = T)
cat("Learning error models is completed.")

# == Check convergence ==
cat("\n")
cat("Forward read error rates (must be decreasing): ", dada2:::checkConvergence(errF), "\n")
cat("Reverse read error rates (must be decreasing): ", dada2:::checkConvergence(errR), "\n")

cat(sprintf("Exporting error plots to `%s`\n",
            file.path(plots.dir, "error_plots.pdf")))
pdf(file = file.path(plots.dir, "error_plots.pdf"))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()
cat("Done.\n")

# == Sample Inference ==
cat("\n")
cat("Sample inference started (it will take a while).\n")
dadaFs <- dada(derepFs, err = errF, selfConsist = TRUE, pool = TRUE,
               multithread = TRUE, verbose = T)
dadaRs <- dada(derepRs, err = errR, selfConsist = TRUE, pool = TRUE,
               multithread = TRUE, verbose = T)
cat("Sample inference is completed.\n")

# == Merge paired reads ==
cat("\n")
cat("Read merging started.\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                      verbose = TRUE, minOverlap = 12, maxMismatch=1)
cat("Read merging is completed.\n")

# == Construct sequence table ==
seqtab <- makeSequenceTable(mergers)

# == Remove chimeras ==
cat("\n")
cat("Chimera removing started (it will take a while).\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
cat("Chimera removing is completed.\n")

# == Track reads through the pipeline, like in tutorial ==
cat("\n")
cat("Tracking reads through the pipeline:")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(head(track))

# == Assign taxonomy ==
cat("\n\n")
cat("Taxonomic classification started (it will take a while).\n")
taxtab <- assignTaxonomy(seqtab.nochim, refFasta = ref.tax.file,
                         multithread = TRUE, verbose = TRUE)
colnames(taxtab) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
cat("\nTaxonomic classification is completed.\n")

# == Give our seq headers more manageable names (ASV_1, ASV_2...) and save seqs to fasta == 
avs_fa_fpath <- file.path(outfiles.dir, "ASV.fa")
cat("\n")
cat(sprintf("Saving ASV sequences to `%s`\n", avs_fa_fpath))
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
# Save sequences in fasta file
write(asv_fasta, avs_fa_fpath)
rm(avs_fa_fpath)

# == Create and save count table ==
avs_count_fpath <- file.path(outfiles.dir, "ASV_abundance_table.tsv")
cat("\n")
cat(sprintf("Saving ASV abundance table to `%s`\n", avs_count_fpath))
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, avs_count_fpath,
            sep="\t", quote=F, col.names=NA)
rm(avs_count_fpath)

# == Create and save tax table: ==
avs_tax_fpath <- file.path(outfiles.dir, "ASVs_taxonomy.tsv")
cat("\n")
cat(sprintf("Saving taxonomy table to `%s`\n", avs_tax_fpath))
asv_tax <- taxtab
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, avs_tax_fpath,
            sep="\t", quote=F, col.names=NA)
rm(avs_tax_fpath)

# == Create a basic phylogenetic tree ==
cat("\n")
cat("Start making basic phylogenetic tree.\n")
seqs <- getSequences(seqtab.nochim)
names(seqs) <- sub(">", "", asv_headers)
# Perform Multiple sequence alignment (MSA)
cat("Performing multiple sequence alignment (it will take a while).\n")
mult.alignment <- msa(seqs, method = "Muscle", type = "dna",
                      order = "input", verbose = TRUE)
cat("Sequences are aligned.")
cat("Now we will create phylogenetic tree and fit it with maximum likelihood method.\n")
phang.align <- as.phyDat(mult.alignment, type = "DNA",
                         names = sub(">", "", asv_headers))
dist.matr <- dist.ml(phang.align) # distance matrix
treeNJ <- NJ(dist.matr) # make starting tree with Neighbor-Joining algorithm
# Update tree -- fit it with maximum likelihood method
fit <- pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 1))
cat("Basic phylogenetic tree is created.\n")

# == Combine data into a phyloseq object ==

colnames(seqtab.nochim) <- sub(">", "", asv_headers)
rownames(taxtab) <- sub(">", "", asv_headers)

cat("\n")
cat("Creating phyloseq object...\n")
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(metadata.df),
               tax_table(taxtab), phy_tree(fitGTR$tree))
ps_merged <- merge_phyloseq(ps, metadata.df)
cat("Phyloseq object is created.\n")

ps.save.fpath <- file.path(outfiles.dir, "phyloseq_object.rds")
cat(sprintf("\nSaving created phyloseq object to `%s`.\n", ps.save.fpath))
saveRDS(ps_merged, ps.save.fpath)

quit(save = 'no', status = 0)
