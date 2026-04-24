library(readxl)
library(writexl)

# define function
percentile_rank <- function(x, data) {
  sum(data >= x) / length(data) * 100
}


# read in all count tables for K562
K562_WT_unedited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/K562_gDNA_edit_counts.xlsx", sheet = "K562 WT unedited cells")
K562_WT_edited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/K562_gDNA_edit_counts.xlsx", sheet = "K562 WT edited cells")
K562_with_marker_unedited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/K562_gDNA_edit_counts.xlsx", sheet = "K562 with marker unedited cells")
K562_with_marker_edited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/K562_gDNA_edit_counts.xlsx", sheet = "K562 with marker edited cells")
K562_without_marker_unedited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/K562_gDNA_edit_counts.xlsx", sheet = "K562 without marker unedited ce")
K562_without_marker_edited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/K562_gDNA_edit_counts.xlsx", sheet = "K562 without marker edited cell")


K562_WT_unedited_cells <- as.data.frame(K562_WT_unedited_cells[-1, c(1,5,7,9)])
K562_WT_unedited_cells$WT_unedited_cells_mean <- rowMeans(K562_WT_unedited_cells[,c(2:4)])

K562_WT_edited_cells <- as.data.frame(K562_WT_edited_cells[-1, c(1,5,7,9)])
K562_WT_edited_cells$WT_edited_cells_mean <- rowMeans(K562_WT_edited_cells[,c(2:4)])

K562_with_marker_unedited_cells <- as.data.frame(K562_with_marker_unedited_cells[-1, c(1,5,7,9)])
K562_with_marker_unedited_cells$with_marker_unedited_cells_mean <- rowMeans(K562_with_marker_unedited_cells[,c(2:4)])

K562_with_marker_edited_cells <- as.data.frame(K562_with_marker_edited_cells[-1, c(1,5,7,9)])
K562_with_marker_edited_cells$with_marker_edited_cells_mean <- rowMeans(K562_with_marker_edited_cells[,c(2:4)])

K562_without_marker_unedited_cells <- as.data.frame(K562_without_marker_unedited_cells[-1, c(1,5,7,9)])
K562_without_marker_unedited_cells$without_marker_unedited_cells_mean <- rowMeans(K562_without_marker_unedited_cells[,c(2:4)])

K562_without_marker_edited_cells <- as.data.frame(K562_without_marker_edited_cells[-1, c(1,5,7,9)])
K562_without_marker_edited_cells$without_marker_edited_cells_mean <- rowMeans(K562_without_marker_edited_cells[,c(2:4)])

# K562 with marker 
K562_with_marker_and_WT <- cbind(K562_with_marker_unedited_cells[c(1,5)], K562_with_marker_edited_cells[5], K562_WT_unedited_cells[5], K562_WT_edited_cells[5])


K562_with_marker_and_WT_odds_ratio <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(K562_with_marker_and_WT_odds_ratio) <- c("Sample ID", "odds ratio", "pVal", "log10(odds ratio)")
for(i in 1:nrow(K562_with_marker_and_WT)){
  table1 <- matrix(c(round(K562_with_marker_and_WT[i,3])+1, round(K562_with_marker_and_WT[i,5])+1, round(K562_with_marker_and_WT[i,2])+1, round(K562_with_marker_and_WT[i,4]+1)), nrow = 2, byrow = TRUE,
                   dimnames = list(cells = c("edited", "unedited"),
                                   Marker = c("with marker", "WT")))
  fisher_results <- fisher.test(table1)
  
  K562_with_marker_and_WT_odds_ratio[i,1] <- K562_with_marker_and_WT[i,1]
  K562_with_marker_and_WT_odds_ratio[i,2] <- fisher_results$estimate
  K562_with_marker_and_WT_odds_ratio[i,3] <- fisher_results$p.value
  K562_with_marker_and_WT_odds_ratio[i,4] <- log10(fisher_results$estimate)
}

K562_with_marker_and_WT_odds_ratio <- K562_with_marker_and_WT_odds_ratio[order(K562_with_marker_and_WT_odds_ratio$`log10(odds ratio)`),]
K562_with_marker_and_WT_odds_ratio$percentile_rank <- sapply(K562_with_marker_and_WT_odds_ratio$`log10(odds ratio)`, percentile_rank, data = K562_with_marker_and_WT_odds_ratio$`log10(odds ratio)`)

write_xlsx(K562_with_marker_and_WT_odds_ratio, "~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Odds_Ratio_and_gDNA/K562_with_marker_and_WT_odds_ratio.xlsx")


# K562 without marker 
K562_without_marker_and_WT <- cbind(K562_without_marker_unedited_cells[c(1,5)], K562_without_marker_edited_cells[5], K562_WT_unedited_cells[5], K562_WT_edited_cells[5])


K562_without_marker_and_WT_odds_ratio <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(K562_without_marker_and_WT_odds_ratio) <- c("Sample ID", "odds ratio", "pVal", "log10(odds ratio)")
for(i in 1:nrow(K562_without_marker_and_WT)){
  table1 <- matrix(c(round(K562_without_marker_and_WT[i,3])+1, round(K562_without_marker_and_WT[i,5])+1, round(K562_without_marker_and_WT[i,2])+1, round(K562_without_marker_and_WT[i,4]+1)), nrow = 2, byrow = TRUE,
                   dimnames = list(cells = c("edited", "unedited"),
                                   Marker = c("with marker", "WT")))
  fisher_results <- fisher.test(table1)
  
  K562_without_marker_and_WT_odds_ratio[i,1] <- K562_without_marker_and_WT[i,1]
  K562_without_marker_and_WT_odds_ratio[i,2] <- fisher_results$estimate
  K562_without_marker_and_WT_odds_ratio[i,3] <- fisher_results$p.value
  K562_without_marker_and_WT_odds_ratio[i,4] <- log10(fisher_results$estimate)
}

K562_without_marker_and_WT_odds_ratio <- K562_without_marker_and_WT_odds_ratio[order(K562_without_marker_and_WT_odds_ratio$`log10(odds ratio)`),]
K562_without_marker_and_WT_odds_ratio$percentile_rank <- sapply(K562_without_marker_and_WT_odds_ratio$`log10(odds ratio)`, percentile_rank, data = K562_without_marker_and_WT_odds_ratio$`log10(odds ratio)`)

write_xlsx(K562_without_marker_and_WT_odds_ratio, "~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Odds_Ratio_and_gDNA/K562_without_marker_and_WT_odds_ratio.xlsx")



# read in all count tables for KBM7
KBM7_WT_unedited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/KBM7_gDNA_edit_counts.xlsx", sheet = "KBM7 WT unedited cells")
KBM7_WT_edited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/KBM7_gDNA_edit_counts.xlsx", sheet = "KBM7 WT edited cells")
KBM7_with_marker_unedited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/KBM7_gDNA_edit_counts.xlsx", sheet = "KBM7 with marker unedited cells")
KBM7_with_marker_edited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/KBM7_gDNA_edit_counts.xlsx", sheet = "KBM7 with marker edited cells")
KBM7_without_marker_unedited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/KBM7_gDNA_edit_counts.xlsx", sheet = "KBM7 without marker unedited ce")
KBM7_without_marker_edited_cells <- read_xlsx("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Count_tables/KBM7_gDNA_edit_counts.xlsx", sheet = "KBM7 without marker edited cell")


KBM7_WT_unedited_cells <- as.data.frame(KBM7_WT_unedited_cells[-1, c(1,5,7,9)])
KBM7_WT_unedited_cells$WT_unedited_cells_mean <- rowMeans(KBM7_WT_unedited_cells[,c(2:4)])

KBM7_WT_edited_cells <- as.data.frame(KBM7_WT_edited_cells[-1, c(1,5,7,9)])
KBM7_WT_edited_cells$WT_edited_cells_mean <- rowMeans(KBM7_WT_edited_cells[,c(2:4)])

KBM7_with_marker_unedited_cells <- as.data.frame(KBM7_with_marker_unedited_cells[-1, c(1,5,7,9)])
KBM7_with_marker_unedited_cells$with_marker_unedited_cells_mean <- rowMeans(KBM7_with_marker_unedited_cells[,c(2:4)])

KBM7_with_marker_edited_cells <- as.data.frame(KBM7_with_marker_edited_cells[-1, c(1,5,7,9)])
KBM7_with_marker_edited_cells$with_marker_edited_cells_mean <- rowMeans(KBM7_with_marker_edited_cells[,c(2:4)])

KBM7_without_marker_unedited_cells <- as.data.frame(KBM7_without_marker_unedited_cells[-1, c(1,5,7,9)])
KBM7_without_marker_unedited_cells$without_marker_unedited_cells_mean <- rowMeans(KBM7_without_marker_unedited_cells[,c(2:4)])

KBM7_without_marker_edited_cells <- as.data.frame(KBM7_without_marker_edited_cells[-1, c(1,5,7,9)])
KBM7_without_marker_edited_cells$without_marker_edited_cells_mean <- rowMeans(KBM7_without_marker_edited_cells[,c(2:4)])

# KBM7 with marker 
KBM7_with_marker_and_WT <- cbind(KBM7_with_marker_unedited_cells[c(1,5)], KBM7_with_marker_edited_cells[5], KBM7_WT_unedited_cells[5], KBM7_WT_edited_cells[5])

KBM7_with_marker_and_WT_odds_ratio <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(KBM7_with_marker_and_WT_odds_ratio) <- c("Sample ID", "odds ratio", "pVal", "log10(odds ratio)")
for(i in 1:nrow(KBM7_with_marker_and_WT)){
  table1 <- matrix(c(round(KBM7_with_marker_and_WT[i,3])+1, round(KBM7_with_marker_and_WT[i,5])+1, round(KBM7_with_marker_and_WT[i,2])+1, round(KBM7_with_marker_and_WT[i,4]+1)), nrow = 2, byrow = TRUE,
                   dimnames = list(cells = c("edited", "unedited"),
                                   Marker = c("with marker", "WT")))
  fisher_results <- fisher.test(table1)
  
  KBM7_with_marker_and_WT_odds_ratio[i,1] <- KBM7_with_marker_and_WT[i,1]
  KBM7_with_marker_and_WT_odds_ratio[i,2] <- fisher_results$estimate
  KBM7_with_marker_and_WT_odds_ratio[i,3] <- fisher_results$p.value
  KBM7_with_marker_and_WT_odds_ratio[i,4] <- log10(fisher_results$estimate)
}

KBM7_with_marker_and_WT_odds_ratio <- KBM7_with_marker_and_WT_odds_ratio[order(KBM7_with_marker_and_WT_odds_ratio$`log10(odds ratio)`),]
KBM7_with_marker_and_WT_odds_ratio$percentile_rank <- sapply(KBM7_with_marker_and_WT_odds_ratio$`log10(odds ratio)`, percentile_rank, data = KBM7_with_marker_and_WT_odds_ratio$`log10(odds ratio)`)
write_xlsx(KBM7_with_marker_and_WT_odds_ratio, "~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Odds_Ratio_and_gDNA/KBM7_with_marker_and_WT_odds_ratio.xlsx")


# KBM7 without marker 
KBM7_without_marker_and_WT <- cbind(KBM7_without_marker_unedited_cells[c(1,5)], KBM7_without_marker_edited_cells[5], KBM7_WT_unedited_cells[5], KBM7_WT_edited_cells[5])


KBM7_without_marker_and_WT_odds_ratio <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(KBM7_without_marker_and_WT_odds_ratio) <- c("Sample ID", "odds ratio", "pVal", "log10(odds ratio)")
for(i in 1:nrow(KBM7_without_marker_and_WT)){
  table1 <- matrix(c(round(KBM7_without_marker_and_WT[i,3])+1, round(KBM7_without_marker_and_WT[i,5])+1, round(KBM7_without_marker_and_WT[i,2])+1, round(KBM7_without_marker_and_WT[i,4]+1)), nrow = 2, byrow = TRUE,
                   dimnames = list(cells = c("edited", "unedited"),
                                   Marker = c("with marker", "WT")))
  fisher_results <- fisher.test(table1)
  
  KBM7_without_marker_and_WT_odds_ratio[i,1] <- KBM7_without_marker_and_WT[i,1]
  KBM7_without_marker_and_WT_odds_ratio[i,2] <- fisher_results$estimate
  KBM7_without_marker_and_WT_odds_ratio[i,3] <- fisher_results$p.value
  KBM7_without_marker_and_WT_odds_ratio[i,4] <- log10(fisher_results$estimate)
}

KBM7_without_marker_and_WT_odds_ratio <- KBM7_without_marker_and_WT_odds_ratio[order(KBM7_without_marker_and_WT_odds_ratio$`log10(odds ratio)`),]
KBM7_without_marker_and_WT_odds_ratio$percentile_rank <- sapply(KBM7_without_marker_and_WT_odds_ratio$`log10(odds ratio)`, percentile_rank, data = KBM7_without_marker_and_WT_odds_ratio$`log10(odds ratio)`)
write_xlsx(KBM7_without_marker_and_WT_odds_ratio, "~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Fastq/Odds_Ratio_and_gDNA/KBM7_without_marker_and_WT_odds_ratio.xlsx")


