library("tidyverse")
library("limma")
library(DESeq2)
library("RColorBrewer")

#setwd("/mnt/sharc/fastdata/bo1cv/covid_mirna/unmapped_stringent//star.dir/")
setwd("./star.dir/")
data_names_files = list.files(pattern=".Log.final.out$")
data_names <- str_remove_all(data_names_files, pattern = "trimmed-")
data_names <- str_remove(data_names, pattern = ".Log.final.out")
data_samples <- lapply(data_names_files, read.table, sep = "\t", skip = 3,fill = T)
names(data_samples) <- data_names

resume_star_stats <- data_samples %>%
  purrr::reduce(full_join, by = "V1")

colnames(resume_star_stats)[2:ncol(resume_star_stats)] <- data_names

write.table(resume_star_stats, "resume_star_stats.txt", quote = F, 
            row.names = TRUE, sep = "\t")
#resume_star_stats <- read.table("star.dir/resume_star_stats.txt", sep = "\t")

#Mapping unique
resume_star_stats %>% filter(str_detect(V1, " Uniquely mapped reads %")) %>%
  t() -> unique_map
unique_map <- as.data.frame(unique_map[2:nrow(unique_map),])
colnames(unique_map) <- "Percent Uniquely mapped reads"
unique_map <- unique_map %>% 
  mutate(`Percent Uniquely mapped reads` = as.numeric(str_remove(`Percent Uniquely mapped reads`, "%"))) %>%
  rownames_to_column("Sample")
ggplot(unique_map, aes(x = Sample, y = `Percent Uniquely mapped reads`)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=`Percent Uniquely mapped reads`), vjust=1.6, color="white", size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
ggsave("Percent_Uniquely_mapped_reads.png")

#Mapping multiple 
resume_star_stats %>% filter(str_detect(V1, " % of reads mapped to multiple loci")) %>%
  t() -> multi_map
multi_map <- as.data.frame(multi_map[2:nrow(multi_map),])
colnames(multi_map) <- "Percent Multi mapped reads"
multi_map <- multi_map %>% 
  mutate(`Percent Multi mapped reads` = as.numeric(str_remove(`Percent Multi mapped reads`, "%"))) %>%
  rownames_to_column("Sample")
ggplot(multi_map, aes(x = Sample, y = `Percent Multi mapped reads`)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=`Percent Multi mapped reads`), vjust=1.6, color="white", size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
ggsave("Percent_Multi_mapped_reads.png")

#" % of reads mapped to too many loci"
resume_star_stats %>% filter(str_detect(V1, " % of reads mapped to too many loci")) %>%
  t() -> multi_map
multi_map <- as.data.frame(multi_map[2:nrow(multi_map),])
colnames(multi_map) <- "Percent mapped to too many loci"
multi_map <- multi_map %>% 
  mutate(`Percent mapped to too many loci` = as.numeric(str_remove(`Percent mapped to too many loci`, "%"))) %>%
  rownames_to_column("Sample")
ggplot(multi_map, aes(x = Sample, y = `Percent mapped to too many loci`)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=`Percent mapped to too many loci`), vjust=1.6, color="white", size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
ggsave("Percent_mapped_to_too_many_loci.png")

#" % of reads too short"
resume_star_stats %>% filter(str_detect(V1, " % of reads unmapped: too short")) %>%
  t() -> multi_map
multi_map <- as.data.frame(multi_map[2:nrow(multi_map),])
colnames(multi_map) <- "Percent reads too short"
multi_map <- multi_map %>% 
  mutate(`Percent reads too short` = as.numeric(str_remove(`Percent reads too short`, "%"))) %>%
  rownames_to_column("Sample")
ggplot(multi_map, aes(x = Sample, y = `Percent reads too short`)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=`Percent reads too short`), vjust=1.6, color="white", size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
ggsave("Percent_reads_too_short.png")

#Too many mismatches
resume_star_stats %>% filter(str_detect(V1, " % of reads unmapped: too many mismatches")) %>%
  t() -> mis_map
mis_map <- as.data.frame(mis_map[2:nrow(mis_map),])
colnames(mis_map) <- "Percent reads too many mismatches (>1)"
mis_map <- mis_map %>% 
  mutate(`Percent reads too many mismatches (>1)` = as.numeric(str_remove(`Percent reads too many mismatches (>1)`, "%"))) %>%
  rownames_to_column("Sample")
ggplot(mis_map, aes(x = Sample, y = `Percent reads too many mismatches (>1)`)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=`Percent reads too many mismatches (>1)`), vjust=1.6, color="white", size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
ggsave("Percent_too_many_mismatches.png")

resume_star_stats %>% filter(str_detect(V1, " % of reads unmapped: other")) %>%
  t() -> unmapp_other
unmapp_other <- as.data.frame(unmapp_other[2:nrow(unmapp_other),])
colnames(unmapp_other) <- "Percent reads unmapped other"
unmapp_other <- unmapp_other %>% 
  mutate(`Percent reads unmapped other` = as.numeric(str_remove(`Percent reads unmapped other`, "%"))) %>%
  rownames_to_column("Sample")
ggplot(unmapp_other, aes(x = Sample, y = `Percent reads unmapped other`)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=`Percent reads unmapped other`), vjust=1.6, color="white", size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()
ggsave("Percent_unmapped_other.png")
