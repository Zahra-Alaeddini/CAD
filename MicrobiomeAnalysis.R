library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)  # For VST
library(limma)   # For removeBatchEffect
library(DECIPHER)
library(phangorn)
library(Biostrings)
library(ape)

batch_paths <- c("sequence_data/dataset1", "sequence_data/dataset2",
                 "sequence_data/dataset3", "sequence_data/dataset4")

metadata_file <- "sequence_data/metadata.txt"

meta <- read.table(metadata_file, header = TRUE, sep = "\t", row.names = 1)
meta$batch <- as.factor(meta$batch)
meta$group_type <- as.factor(meta$group_type)

seqtabs <- list()
track_list <- list()

for (i in 1:length(batch_paths)) {
  path <- batch_paths[i]
  
  fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))
  
  sample.names <- sapply(strsplit(basename(fnFs), "_R1"), [, 1)
  
  pdf(paste0("quality_profiles_forward_batch", i, ".pdf"))
  
  plotQualityProfile(fnFs[sample(1:length(fnFs), min(10, length(fnFs)))])
  dev.off()
  
  pdf(paste0("quality_profiles_reverse_batch", i, ".pdf"))
  plotQualityProfile(fnRs[sample(1:length(fnRs), min(10, length(fnRs)))])
  dev.off()
  
  filt_path <- file.path(path, "filtered")
  
  if (!dir.exists(filt_path)) dir.create(filt_path)
  filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       maxN = 0, maxEE = c(2, 2), truncQ = 2,
                       minLen = c(250, 220), trimLeft = c(20, 20),
                       rm.phix = TRUE, compress = TRUE, multithread = TRUE)
  
  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)
  
  dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
  
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
  
  seqtab <- makeSequenceTable(mergers)
  
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
  
  seqtabs[[i]] <- seqtab.nochim
  
  getN <- function(x) sum(getUniques(x))
  
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  
  rownames(track) <- sample.names
  
  track_list[[i]] <- track
}

seqtab.all <- mergeSequenceTables(tables = seqtabs)

dim(seqtab.all)

table(nchar(getSequences(seqtab.all)))

#Download from: https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz
silva_ref <- "path/to/silva_nr99_v138.1_train_set.fa.gz"

taxa <- assignTaxonomy(seqtab.all, silva_ref, multithread = TRUE)

seqs <- colnames(seqtab.all)

seqs <- DNAStringSet(seqs)

names(seqs) <- paste0("ASV_", 1:length(seqs))

alignment <- AlignSeqs(seqs, processors = 4)

phydat <- phyDat(as(alignment, "matrix"), type = "DNA")

dm <- dist.ml(phydat)

treeNJ <- NJ(dm)

fit <- pml(treeNJ, data = phydat)

fitGTR <- update(fit, k = 4, inv = 0.2)

fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

tree <- midpoint(fitGTR$tree)

ps <- phyloseq(otu_table(seqtab.all, taxa_are_rows = FALSE),
               tax_table(taxa),
               phy_tree(tree),
               sample_data(meta))

ps <- subset_taxa(ps, Kingdom == "Bacteria")

ps <- subset_taxa(ps, Family != "Mitochondria" & Order != "Chloroplast")

ps <- prune_taxa(taxa_sums(ps) > 1, ps)

ps <- prune_samples(sample_sums(ps) >= 5000, ps)

sample_data(ps) <- sample_data(sample_data(ps)[sample_names(ps), ])

rel_abund <- taxa_sums(ps) / sum(taxa_sums(ps)) * 100

keep_taxa <- rel_abund >= 0.005

ps <- prune_taxa(keep_taxa, ps)

dds <- DESeqDataSetFromPhyloseq(ps)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

uncorrected_counts <- assay(vsd)

dist_uncorrected <- vegdist(uncorrected_counts, method = "bray")

pcoa_unc <- cmdscale(dist_uncorrected, k = 2)

pcoa_unc_df <- data.frame(PC1 = pcoa_unc[,1], PC2 = pcoa_unc[,2], batch = colData(vsd)$batch, group_type = colData(vsd)$group_type)

pdf("beta_bray_pcoa_uncorrected_batch.pdf")

ggplot(pcoa_unc_df, aes(PC1, PC2, color = batch)) + geom_point() + theme_minimal() + ggtitle("PCoA Before Batch Correction - Colored by Batch")
dev.off()

pdf("beta_bray_pcoa_uncorrected_group.pdf")
ggplot(pcoa_unc_df, aes(PC1, PC2, color = group_type)) + geom_point() + theme_minimal() + ggtitle("PCoA Before Batch Correction - Colored by Group Type")
dev.off()

adonis_batch <- adonis2(dist_uncorrected ~ batch, data = as.data.frame(colData(vsd)), permutations = 999)

print("PERMANOVA for Batch Effect Before Correction:")

print(adonis_batch)

adonis_group_unc <- adonis2(dist_uncorrected ~ group_type, data = as.data.frame(colData(vsd)), permutations = 999)

print("PERMANOVA for Group Type Before Correction:")

print(adonis_group_unc)

design <- model.matrix(~ group_type, data = colData(vsd))

corrected_vst <- removeBatchEffect(assay(vsd), batch = vsd$batch, design = design)

assay(vsd) <- corrected_vst

count_matrix <- assay(vsd)

pdf("taxa_bar_plots.pdf")

plot_bar(ps, fill = "Genus")
dev.off()

shannon <- diversity(count_matrix, index = "shannon")

alpha_df <- data.frame(Shannon = shannon, group_type = colData(vsd)$group_type)

pdf("alpha_shannon_significance.pdf")

ggplot(alpha_df, aes(group_type, Shannon)) + geom_boxplot() + theme_minimal()
dev.off()

wilcox.test(Shannon ~ group_type, data = alpha_df)

dist_matrix <- vegdist(count_matrix, method = "bray")

pcoa <- cmdscale(dist_matrix, k = 2)

pcoa_df <- data.frame(PC1 = pcoa[,1], PC2 = pcoa[,2], group_type = colData(vsd)$group_type)

pdf("beta_bray_pcoa.pdf")
ggplot(pcoa_df, aes(PC1, PC2, color = group_type)) + geom_point() + theme_minimal()
dev.off()

adonis2(dist_matrix ~ group_type, data = as.data.frame(colData(vsd)), permutations = 999)