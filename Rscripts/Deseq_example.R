
#Filtering
mirbqse_results <- mirbqse_results[rowSums(mirbqse_results[2:ncol(mirbqse_results)]) >= 1,]
rownames(mirbqse_results) <- mirbqse_results$V1
mirbqse_results$V1 <- NULL
#mirbqse_results <- column_to_rownames(mirbqse_results, "V1")
mirbqse_results <- as.matrix(mirbqse_results)

# write.table(mirbqse_results,"/mnt/sharc/fastdata/bo1cv/covid_mirna/bowtie_mirbase/bowtie.dir/resume_idxstats.tsv",
#             quote = F, sep = "\t", row.names = TRUE)

samples <- data.frame(samples = colnames(mirbqse_results), 
                      condition = str_extract(colnames(mirbqse_results), regex("H|C|A")))

design <- model.matrix(~0+samples$condition)
colnames(design) <- c("A","C","H")

samples$condition <- factor(samples$condition)

dds <- DESeqDataSetFromMatrix(mirbqse_results,
                              colData = samples,
                              design = ~ condition)

dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE, nsub = sum( rowMeans( counts(dds, normalized=TRUE)) > 5 )) 

pcaData <- plotPCA(vsd, intgroup=c("condition", "samples"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#jpeg("plotPCA.png")
ggplot(pcaData, aes(PC1, PC2, color=samples, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#dev.off()


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$samples, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#jpeg("plotDist.png")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)