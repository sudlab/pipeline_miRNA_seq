library("tidyverse")
library("limma")
library("DESeq2")
library("RColorBrewer")
library("optparse")
#Options parser
option_list = list(
  make_option(c("-o", "--output"),
              type="character",
              dest = "output",
              help="output name")
)
arguments <- parse_args(OptionParser(option_list = option_list))

# arguments <- data.frame(output = "counts_first.dir/Primary_resume_counts.txt")
# setwd("/mnt/sharc/shared/sudlab1/General/projects/Favour_EKT/pipeline_seq_favour/counts_first.dir/")

the_pattern <- paste0("counts",str_extract(arguments$output, "MultiFrac|Multi|Primary|Overlap|Unique"),".txt$")
data_names_files = list.files(pattern=the_pattern, recursive = TRUE)
data_names <- str_remove_all(data_names_files, pattern = "trimmed-")
data_names <- str_remove(data_names, pattern = "_counts[:alpha:]+.txt")
data_names <- str_remove(data_names, pattern = "counts.+.dir/")
data_samples <- lapply(data_names_files, read.table, header = T)
names(data_samples) <- data_names

resume_data <- data_samples %>%
purrr::reduce(full_join, by = c("Geneid","Chr","Start","End","Length","Strand"))

colnames(resume_data)[7:ncol(resume_data)] <- data_names

resume_data <- resume_data %>%
  group_by(Geneid) %>%
  summarise_at(colnames(resume_data)[7:ncol(resume_data)] , sum)

#Due to GTF from miRBase
if (str_detect(resume_data$Geneid[1], "transcript")) {
  resume_data <- resume_data %>%
    mutate(Geneid = str_remove_all(Geneid,"transcript:"))
}

#Filter at least one read
resume_data <- resume_data[rowSums(resume_data[2:ncol(resume_data)]) >= 1,]
resume_data <- column_to_rownames(resume_data, "Geneid")
resume_data <- as.matrix(resume_data)


##PCA an dsample dist with deseq
#Sample table for deseq
samples <- data.frame(samples = colnames(resume_data), 
                      condition = str_extract(colnames(resume_data), regex("OE|CTRL")),
                      day = str_extract(colnames(resume_data), regex("D.")),
                      replicate = str_extract(colnames(resume_data), regex("..Unmapped|.$")))
samples$condition <- factor(samples$condition)

dds <- try(DESeqDataSetFromMatrix(resume_data,
                              colData = samples,
                              design = ~ condition))

if (class(dds) == "try-error") {
  dds <- try(DESeqDataSetFromMatrix(round(resume_data),
                                    colData = samples,
                                    design = ~ condition))
}

dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE, nsub = sum( rowMeans( counts(dds, normalized=TRUE)) > 5 )) 

pcaData <- plotPCA(vsd, intgroup=c("condition", "samples"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

jpeg(paste0("plotPCA",str_extract(arguments$output, "MultiFrac|Multi|Primary|Overlap"),".png"))
ggplot(pcaData, aes(PC1, PC2, color=samples, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()


# sampleDists <- dist(t(assay(vsd)))
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$samples, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# jpeg(paste0("plotDist",str_extract(arguments$output, "MultiFrac|Multi|Primary|Overlap"),".png"))
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
# dev.off()

##

write.table(resume_data, arguments$output, quote = F, 
              row.names = TRUE, sep = "\t")

