rm(list=ls())
library(readxl)
library(writexl)
library(stringi)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)


setwd("~/Design_pegRNAs_additional_STOP_codon")

# Input sequence of the gene 
gene_seq <- read_file("Chr17_full_sequence.txt")
nchar(gene_seq)


# read in coordinates of exons 
exons <- as.data.frame(read.table("NF1-Exons-positions.txt"))
exons <- exons[-c(1,2),]
colnames(exons) <- c("Exon", "Start", "End")
exons <- as.data.frame(sapply(exons,as.numeric))
exon27 <- filter(exons, exons$Exon == 27)
exon27_start <- exon27[1,2]
exon27_end <- exon27[1,3]
exon27_and_intron_seq <- str_sub(gene_seq,exon27_start - 360, exon27_end + 360)


# function to do synonymous mutations
get_synonym_mutation <- function(codon){
  if(codon == "ATT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "ATC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "ATA"){
    ref <- "A"
    var <- "T"
    return(var) 
  }
  else if(codon == "ATG"){
    as <- "ATG"
    return(as)
  }
  else if(codon == "ACT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "ACC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "ACA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "ACG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "AAT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "AAC"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "AAA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "AAG"){
    ref <- "G"
    var <- "A"
    return(var) 
  }
  else if(codon == "AGT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "AGC"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "AGA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "AGG"){
    ref <- "G"
    var <- "A"
    return(var) 
  }
  else if(codon == "GTT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "GTC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "GTA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "GTG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "GCT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "GCC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "GCA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "GCG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "GAT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "GAC"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "GAA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "GAG"){
    ref <- "G"
    var <- "A"
    return(var) 
  }
  else if(codon == "GGT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "GGC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "GGA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "GGG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "TTT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "TTC"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "TTA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "TTG"){
    ref <- "G"
    var <- "A"
    return(var) 
  }
  else if(codon == "TCT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "TCC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "TCA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "TCG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "TAT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "TAC"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "TAA"){
    as <- "TAA"
    return(as)
  }
  else if(codon == "TAG"){
    as <- "TAG"
    return(as)
  }
  else if(codon == "TGT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "TGC"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "TGG"){
    as <- "TGG"
    return(as)
  }
  else if(codon == "TGA"){
    as <- "TGA"
    return(as)
  }
  else if(codon == "CTT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "CTC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "CTA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "CTG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "CCT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "CCC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "CCA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "CCG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "CAT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "CAC"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "CAA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "CAG"){
    ref <- "G"
    var <- "A"
    return(var) 
  }
  else if(codon == "CGT"){
    ref <- "T"
    var <- "C"
    return(var) 
  }
  else if(codon == "CGC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "CGA"){
    ref <- "A"
    var <- "G"
    return(var) 
  }
  else if(codon == "CGG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
}

# function to get stop codon
get_stop_codon_pos1 <- function(codon){
  ### on position 1 of a codon
  if(codon == "CAA"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "AAA"){
    ref <- "A"
    var <- "T"
    return(var) 
  }
  else if(codon == "GAA"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "CAG"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "AAG"){
    ref <- "A"
    var <- "T"
    return(var) 
  }
  else if(codon == "GAG"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
  else if(codon == "CGA"){
    ref <- "C"
    var <- "T"
    return(var) 
  }
  else if(codon == "AGA"){
    ref <- "A"
    var <- "T"
    return(var) 
  }
  else if(codon == "GGA"){
    ref <- "G"
    var <- "T"
    return(var) 
  }
}
check_pos1 <- c("CAA", "AAA", "GAA", "CAG", "AAG", "GAG", "CGA", "AGA", "GGA")

get_stop_codon_pos2 <- function(codon){
  ### on position 2 of a codon
  if(codon == "TCA"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "TTA"){
    ref <- "T"
    var <- "A"
    return(var) 
  }
  else if(codon == "TGG"){
    ref <- "G"
    var <- "A"
    return(var) 
  }
  else if(codon == "TCG"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "TTG"){
    ref <- "T"
    var <- "A"
    return(var) 
  }
}
check_pos2 <- c("TCA", "TTA", "TGG", "TCG", "TTG")

get_stop_codon_pos3 <- function(codon){
  ### on position 3 of a codon
  if(codon == "TAC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
  else if(codon == "TAT"){
    ref <- "T"
    var <- "A"
    return(var) 
  }
  
  else if(codon == "TGT"){
    ref <- "T"
    var <- "A"
    return(var) 
  }
  else if(codon == "TGC"){
    ref <- "C"
    var <- "A"
    return(var) 
  }
}
check_pos3 <- c("TAC", "TAT", "TGT", "TGC")

##################################################################################################################

## create input list for additional nonsense mutations
input_list <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(input_list) <- c("InputID","Codon_pos", "Codon_count", "MutationType", "InputSequence") 
exon_edit_subseq <- exon27_and_intron_seq
t = 0
counter = 1
codon_end = 362
codon_counter = 1
for(i in 363:nchar(exon27_and_intron_seq)){ # starts two nucleotides after the beginning of exon 27 as exon 27 starts with an incomplete codon
  exon_edit_subseq <- exon27_and_intron_seq
  if(i >= 573){
    break
  }
  else if(i <= 572 && counter == 3){
    counter = 1
    codon_counter = codon_counter + 1
    codon_end = codon_end + 3
    nt_1 <- codon_end - 2
    nt_2 <- codon_end - 1
    nt_3 <- codon_end 
    codon <- str_sub(exon27_and_intron_seq, nt_1, nt_3)
    if(codon %in% check_pos1){
      var <- get_stop_codon_pos1(codon)
      edit_pos <- nt_1 - 363 + 3 # nt_1 because edit is on first position; -363 + 3 to get position on exon 27 sequence only
      edit_pos_genome <- 29560020 + edit_pos -1
      ref <- str_sub(codon, 1, 1)
      edit <- str_c(var) 
      count_char_edit <- nchar(edit)
      str_sub(exon_edit_subseq, nt_1, nt_1) <- edit
      codon_before <- nt_1 - 3
      codon_after <- nt_3 + 3
      if(codon_after <= 574){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 + 3, nt_3 + 3)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_1 + count_char_edit + 4, nt_1 + count_char_edit + 4) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - 10, nt_1 + count_char_edit + 4 + count_char_edit_synonym -1 + 10)
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- exon_edit_subseq
          
        }
        else if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          if(codon_before >= 363){
            codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 - 3, nt_1 - 1)
            if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
              var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
              ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
              edit_synonym <- str_c(var_synonym)
              count_char_edit_synonym <- nchar(edit_synonym)
              str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
              exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_1 + count_char_edit + 10)
              t = t+1
              input_list[t,1] <- paste0("STOP_", edit_pos_genome)
              input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list[t,3] <- codon_counter
              input_list[t,4] <- "Nonsense mutation"
              input_list[t,5] <- exon_edit_subseq
            }
            if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
              t = t+1
              input_list[t,1] <- paste0("STOP_", edit_pos_genome)
              input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list[t,3] <- codon_counter
              input_list[t,4] <- "Nonsense mutation"
              input_list[t,5] <- "Codon before and after edit is ATG, TGG, TGA, TAA or TAG."
            }
          }
          else if(codon_before < 363){
            t = t+1
            input_list[t,1] <- paste0("STOP_", edit_pos_genome)
            input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
            input_list[t,3] <- codon_counter
            input_list[t,4] <- "Nonsense mutation"
            input_list[t,5] <- "Codon after is ATG, TGG, TGA, TAA or TAG and there is no codon before."
          }
        }
      }
      else if(codon_after > 574 & codon_before >= 363){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 -3, nt_1 -1)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_1 + count_char_edit + 10)
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- exon_edit_subseq
        }
        if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- "Codon before is ATG, TGG, TGA, TAA or TAG and there is no codon after."
        }
      }
    }
    else if(codon %in% check_pos2){
      var <- get_stop_codon_pos2(codon)
      edit_pos <- nt_2 - 363 + 3 # nt_2 because edit is on second position; -363 + 3 to get postion on exon 27 sequence only
      edit_pos_genome <- 29560020 + edit_pos -1
      ref <- str_sub(codon, 2, 2)
      edit <- str_c(var) 
      count_char_edit <- nchar(edit)
      str_sub(exon_edit_subseq, nt_2, nt_2) <- edit
      codon_before <- nt_1 - 3
      codon_after <- nt_3 + 3
      if(codon_after <= 574){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 + 3, nt_3 + 3)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_2 + count_char_edit + 3, nt_2 + count_char_edit + 3) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_2 - 10, nt_2 + count_char_edit + 3 + count_char_edit_synonym -1 + 10)
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- exon_edit_subseq
        }
        else if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          if(codon_before >= 363){
            codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 -3, nt_1 -1)
            if(codon_for_synonym_mut != "ATG" | codon_for_synonym_mut != "TGG" | codon_for_synonym_mut != "TGA" | codon_for_synonym_mut != "TAA" | codon_for_synonym_mut != "TAG"){
              var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
              ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
              edit_synonym <- str_c(var_synonym)
              count_char_edit_synonym <- nchar(edit_synonym)
              str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
              exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_2 + count_char_edit + 10)
              t = t+1
              input_list[t,1] <- paste0("STOP_", edit_pos_genome)
              input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list[t,3] <- codon_counter
              input_list[t,4] <- "Nonsense mutation"
              input_list[t,5] <- exon_edit_subseq
            }
            if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
              t = t+1
              input_list[t,1] <- paste0("STOP_", edit_pos_genome)
              input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list[t,3] <- codon_counter
              input_list[t,4] <- "Nonsense mutation"
              input_list[t,5] <- "Codon before and after edit is ATG, TGG, TGA, TAA or TAG."
            }
          }
          else if(codon_before < 363){
            t = t+1
            input_list[t,1] <- paste0("STOP_", edit_pos_genome)
            input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
            input_list[t,3] <- codon_counter
            input_list[t,4] <- "Nonsense mutation"
            input_list[t,5] <- "Codon after is ATG, TGG, TGA, TAA or TAG and there is no codon before."
          }
        }
      }
      else if(codon_after > 574 & codon_before >= 363){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 -3, nt_1 -1)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_2 + count_char_edit + 10)
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- exon_edit_subseq
        }
        if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- "Codon before is ATG, TGG, TGA, TAA or TAG and there is no codon after."
        }
      }
    }
    else if(codon %in% check_pos3){
      var <- get_stop_codon_pos3(codon)
      edit_pos <- nt_3 - 363 + 3 # nt_3 because edit is on third position; -363 + 3 to get postion on exon 27 sequence only
      edit_pos_genome <- 29560020 + edit_pos -1
      ref <- str_sub(codon, 3, 3)
      edit <- str_c(var) 
      count_char_edit <- nchar(edit)
      str_sub(exon_edit_subseq, nt_3, nt_3) <- edit
      codon_before <- nt_1 - 3
      codon_after <- nt_3 + 3
      if(codon_after <= 574){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 + 3, nt_3 + 3)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_3 + count_char_edit + 2, nt_3 + count_char_edit + 2) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_3 - 10, nt_3 + count_char_edit + 2 + count_char_edit_synonym -1 + 10)
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- exon_edit_subseq
        }
        else if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          if(codon_before >= 363){
            codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 -3, nt_1 -1)
            if(codon_for_synonym_mut != "ATG" | codon_for_synonym_mut != "TGG" | codon_for_synonym_mut != "TGA" | codon_for_synonym_mut != "TAA" | codon_for_synonym_mut != "TAG"){
              var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
              ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
              edit_synonym <- str_c(var_synonym)
              count_char_edit_synonym <- nchar(edit_synonym)
              str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
              exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_2 + count_char_edit + 10)
              t = t+1
              input_list[t,1] <- paste0("STOP_", edit_pos_genome)
              input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list[t,3] <- codon_counter
              input_list[t,4] <- "Nonsense mutation"
              input_list[t,5] <- exon_edit_subseq
            }
            if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
              t = t+1
              input_list[t,1] <- paste0("STOP_", edit_pos_genome)
              input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list[t,3] <- codon_counter
              input_list[t,4] <- "Nonsense mutation"
              input_list[t,5] <- "Codon before and after edit is ATG, TGG, TGA, TAA or TAG."
            }
          }
          else if(codon_before < 363){
            t = t+1
            input_list[t,1] <- paste0("STOP_", edit_pos_genome)
            input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
            input_list[t,3] <- codon_counter
            input_list[t,4] <- "Nonsense mutation"
            input_list[t,5] <- "Codon after is ATG, TGG, TGA, TAA or TAG and there is no codon before."
          }
        }
      }
      else if(codon_after > 574 & codon_before >= 363){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 -3, nt_1 -1)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_3 + count_char_edit + 10)
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- exon_edit_subseq
        }
        if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          t = t+1
          input_list[t,1] <- paste0("STOP_", edit_pos_genome)
          input_list[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list[t,3] <- codon_counter
          input_list[t,4] <- "Nonsense mutation"
          input_list[t,5] <- "Codon before is ATG, TGG, TGA, TAA or TAG and there is no codon after."
        }
      }
    }
    else{
      counter = 1
      i = i+1
    }
  }
  else if(i <= 572 && counter < 3){
    counter = counter + 1
    i = i+1
  }
}


#write_xlsx(input_list, "Searching_sequences_stop_gDNA.xlsx") 


### design searching sequences for additional synonymous mutations

input_list_control <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(input_list_control) <- c("InputID","Codon_pos", "Codon_count", "MutationType", "InputSequence") 
exon_edit_subseq <- exon27_and_intron_seq
t = 0
counter = 1
codon_end = 362
codon_counter = 1
for(i in 363:nchar(exon27_and_intron_seq)){ # starts two nucleotides after the beginning of exon 27 as exon 27 starts with an incomplete exon
  exon_edit_subseq <- exon27_and_intron_seq
  if(i >= 573){
    break
  }
  else if(i <= 572 && counter == 3){
    counter = 1
    codon_counter = codon_counter + 1
    codon_end = codon_end + 3
    nt_1 <- codon_end - 2
    nt_2 <- codon_end - 1
    nt_3 <- codon_end 
    codon <- str_sub(exon27_and_intron_seq, nt_1, nt_3)
    var <- get_synonym_mutation(codon)
    edit_pos <- nt_3 - 363 + 3
    edit_pos_genome <- 29560020 + edit_pos -1
    ref <- str_sub(codon, 3, 3)
    edit <- str_c(var) 
    count_char_edit <- nchar(edit)
    if(var != "ATG" & var != "TGG" & var != "TGA" & var != "TAA" & var != "TAG"){
      str_sub(exon_edit_subseq, nt_3, nt_3) <- edit
      codon_before <- nt_1 - 3
      codon_after <- nt_3 + 3
      if(codon_after <= 574){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 + 3, nt_3 + 3)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_3 + count_char_edit + 2, nt_3 + count_char_edit + 2) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_3 - 10, nt_3 + count_char_edit + 2 + count_char_edit_synonym -1 + 10)
          t = t+1
          input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos_genome)
          input_list_control[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list_control[t,3] <- codon_counter
          input_list_control[t,4] <- "Synonymous mutation"
          input_list_control[t,5] <- exon_edit_subseq
        }
        else if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          if(codon_before >= 363){
            codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 -3, nt_1 -1)
            if(codon_for_synonym_mut != "ATG" | codon_for_synonym_mut != "TGG" | codon_for_synonym_mut != "TGA" | codon_for_synonym_mut != "TAA" | codon_for_synonym_mut != "TAG"){
              var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
              ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
              edit_synonym <- str_c(var_synonym)
              count_char_edit_synonym <- nchar(edit_synonym)
              str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
              exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_2 + count_char_edit + 10)
              t = t+1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos_genome)
              input_list_control[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list_control[t,3] <- codon_counter
              input_list_control[t,4] <- "Synonymous mutation"
              input_list_control[t,5] <- exon_edit_subseq
            }
            if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
              t = t+1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos_genome)
              input_list_control[t,2] <- paste0("nucleotide_pos_", edit_pos)
              input_list_control[t,3] <- codon_counter
              input_list_control[t,4] <- "Synonymous mutation"
              input_list_control[t,5] <- "Codon before and after edit is ATG, TGG, TGA, TAA or TAG."
            }
          }
          else if(codon_before < 363){
            t = t+1
            input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos_genome)
            input_list_control[t,2] <- paste0("nucleotide_pos_", edit_pos)
            input_list_control[t,3] <- codon_counter
            input_list_control[t,4] <- "Synonymous mutation"
            input_list_control[t,5] <- "Codon after is ATG, TGG, TGA, TAA or TAG and there is no codon before."
          }
        }
      }
      else if(codon_after > 574 & codon_before >= 363){
        codon_for_synonym_mut <- str_sub(exon27_and_intron_seq, nt_1 -3, nt_1 -1)
        if(codon_for_synonym_mut != "ATG" & codon_for_synonym_mut != "TGG" & codon_for_synonym_mut != "TGA" & codon_for_synonym_mut != "TAA" & codon_for_synonym_mut != "TAG"){
          var_synonym <- get_synonym_mutation(codon_for_synonym_mut)
          ref_synonym <- str_sub(codon_for_synonym_mut,3,3)
          edit_synonym <- str_c(var_synonym)
          count_char_edit_synonym <- nchar(edit_synonym)
          str_sub(exon_edit_subseq, nt_1 - 1, nt_1 -1) <- edit_synonym
          exon_edit_subseq <- str_sub(exon_edit_subseq, nt_1 - count_char_edit_synonym - 10, nt_3 + count_char_edit + 10)
          t = t+1
          input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos_genome)
          input_list_control[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list_control[t,3] <- codon_counter
          input_list_control[t,4] <- "Synonymous mutation"
          input_list_control[t,5] <- exon_edit_subseq
        }
        if(codon_for_synonym_mut == "ATG" | codon_for_synonym_mut == "TGG" | codon_for_synonym_mut == "TGA" | codon_for_synonym_mut == "TAA" | codon_for_synonym_mut == "TAG"){
          t = t+1
          input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos_genome)
          input_list_control[t,2] <- paste0("nucleotide_pos_", edit_pos)
          input_list_control[t,3] <- codon_counter
          input_list_control[t,4] <- "Synonymous mutation"
          input_list_control[t,5] <- "Codon before is ATG, TGG, TGA, TAA or TAG and there is no codon after."
        }
      }
      else{
        counter = 1
        i = i+1
      }
    }
    else if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
      counter = 1
      i = i+1
    }
  }
  else if(i <= 572 && counter < 3){
    counter = counter + 1
    i = i+1
  }
}


#write_xlsx(input_list_control, "Searching_sequences_synonymous_edit_stop_gDNA.xlsx") 

input_searching <- rbind(input_list[,c(1,3,5)], input_list_control[,c(1,3,5)])
write_xlsx(input_searching, "Searching_sequences_stop_and_synonymous_edit_gDNA.xlsx")
