rm(list=ls())
library(readxl)
library(writexl)
library(stringi)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)


setwd("~/Design_pegRNAs_cBioPortal_mutations")

# Input sequence of the gene 
gene_seq <- read_file("Chr17_full_sequence.txt")
nchar(gene_seq)

# read in coordinates of exons 
exons <- as.data.frame(read.table("NF1-Exons-positions.txt"))
exons <- exons[-c(1,2),]
colnames(exons) <- c("Exon", "Start", "End")
exons <- as.data.frame(sapply(exons,as.numeric))

# function to do synonymous mutations
get_silent_mutation <- function(codon){
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

##################################################################################################################

## create input list for synoynmous mutations
input_list_control <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(input_list_control) <- c("InputID","Codon_pos", "MutationType", "InputSequence") 
exon_edit <- gene_seq
t = 0
codon_counter = 0
exon_seq_2=""
exon_seq_1=""
for (i in 1:nrow(exons)){
  exon_edit <- gene_seq
  exon_start_position <- exons[i,2]
  exon_end_position <- exons[i,3]
  exon_seq <- str_sub(gene_seq, exon_start_position, exon_end_position)
  length_exon <- nchar(exon_seq)
  if(exon_seq_1 != exon_seq && exon_seq_2 != exon_seq){
    counter = 1
    codon_end = 0
    for (j in 1:nchar(exon_seq)){
      exon_edit_subseq <- exon_seq
      marker_codon_seq <- gene_seq
      exon_edit <- gene_seq
      if(length_exon >= j && counter == 3){
        counter = 1
        codon_end = codon_end + 3
        nt_1 <- codon_end - 2
        nt_2 <- codon_end - 1
        nt_3 <- codon_end 
        codon <- str_sub(exon_seq, nt_1, nt_3)
        var <- get_silent_mutation(codon)
        if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
          edit_pos <- exon_start_position  + nt_2
          text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
          t=t+1
          codon_counter = codon_counter + 1
          input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
          input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
          input_list_control[t,3] <- "Synonymous"
          input_list_control[t,4] <- text
        }
        else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
          ref <- str_sub(codon, 3, 3)
          edit <- var 
          str_sub(exon_edit_subseq, nt_3, nt_3) <- edit
          str_sub(exon_edit, exon_start_position, exon_end_position) <- exon_edit_subseq
          edit_pos <- exon_start_position  + nt_2
          right_distance <- exon_end_position - edit_pos
          left_distance <- edit_pos - exon_start_position
          if(right_distance >= 3){
            marker_edit_pos <- edit_pos + 3
            marker_codon_start <- edit_pos + 1
            marker_codon_end <- edit_pos + 3
            marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
            marker_ref <- str_sub(marker_codon, 3, 3)
            if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
              marker_var <- get_silent_mutation(marker_codon)
              marker_edit <- marker_var
              str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
              left_distance_2 <- edit_pos - exon_start_position
              right_distance_2 <- exon_end_position - marker_edit_pos
              if(left_distance_2 >= 10 && right_distance_2 >= 10){
                start_position_long <- edit_pos-10
                end_position_long <- marker_codon_end+10
                exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                codon_counter = codon_counter + 1
                t=t+1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- exon_edit
              }
              else if(left_distance_2 >= 10 && right_distance_2 < 10){
                if(i == nrow(exons)){
                  start_position_long <- edit_pos-10
                  end_position_long <- exon_end_position
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(i < nrow(exons)){
                  exon_start_position_after <- exons[i+1,2]
                  exon_end_position_after <- exons[i+1,3]
                  start_position_long <- edit_pos-10
                  end_position_long <- exon_end_position
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  number_of_nts <- 10 - (exon_end_position - marker_codon_end) 
                  exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                  exon_edit <- paste(exon_edit, exon_after, sep = "")
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
              }
              else if(left_distance_2 < 10 && right_distance_2 >= 10){
                if(i == 1){
                  start_position_long <- exon_start_position
                  end_position_long <- marker_codon_end+10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(i > 1){
                  exon_start_position_before <- exons[i-1,2]
                  exon_end_position_before <- exons[i-1,3]
                  start_position_long <- exon_start_position
                  end_position_long <- marker_codon_end + 10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  number_of_nts <- 10- (edit_pos - exon_start_position)
                  exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                  exon_edit <- paste(exon_before, exon_edit, sep = "")
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
              }
            }
            else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
              if(left_distance >= 5){
                marker_edit_pos <- edit_pos - 3
                marker_codon_start <- edit_pos - 5
                marker_codon_end <- edit_pos - 3
                marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
                marker_ref <- str_sub(marker_codon, 3, 3)
                if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                  marker_var <- get_silent_mutation(marker_codon)
                  marker_edit <- marker_var
                  str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                  left_distance_2 <- marker_edit_pos - exon_start_position
                  right_distance_2 <- exon_end_position - edit_pos
                  if(left_distance_2 >= 10 && right_distance_2 >= 10){
                    start_position_long <- marker_codon_end-10
                    end_position_long <- edit_pos+10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(left_distance_2 >= 10 && right_distance_2 < 10){
                    if(i == nrow(exons)){
                      start_position_long <- marker_codon_end-10
                      end_position_long <- exon_end_position
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                    else if(i < nrow(exons)){
                      exon_start_position_after <- exons[i+1,2]
                      exon_end_position_after <- exons[i+1,3]
                      start_position_long <- marker_codon_end-10
                      end_position_long <- exon_end_position
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      number_of_nts <- 10 - (exon_end_position - edit_pos) 
                      exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                      exon_edit <- paste(exon_edit, exon_after, sep = "")
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                  }
                  else if(left_distance_2 < 10 && right_distance_2 >= 10){
                    if(i == 1){
                      start_position_long <- exon_start_position
                      end_position_long <- edit_pos+10
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                    else if(i > 1){
                      exon_start_position_before <- exons[i-1,2]
                      exon_end_position_before <- exons[i-1,3]
                      start_position_long <- exon_start_position
                      end_position_long <- edit_pos + 10
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      number_of_nts <- 10- (marker_codon_end - exon_start_position)
                      exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                      exon_edit <- paste( exon_before, exon_edit, sep = "")
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                  }
                }
                else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                  edit_pos <- exon_start_position  + nt_2
                  text <- "ATG, TGG, TAA, TAG and TGA are before and after the edit."
                  t=t+1
                  codon_counter = codon_counter + 1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- text
                }
              }
              else if(left_distance < 5){
                edit_pos <- exon_start_position  + nt_2
                text <- "ATG, TGG, TAA, TAG and TGA are after the edit and not enough bp before the edit."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
          }
          else if(right_distance < 3){
            if(left_distance >= 5){
              marker_edit_pos <- edit_pos - 3
              marker_codon_start <- edit_pos - 5
              marker_codon_end <- edit_pos - 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                left_distance_2 <- marker_edit_pos - exon_start_position
                right_distance_2 <- exon_end_position - edit_pos
                if(left_distance_2 >= 10 && right_distance_2 >= 10){
                  start_position_long <- marker_codon_end-10
                  end_position_long <- edit_pos+10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(left_distance_2 >= 10 && right_distance_2 < 10){
                  if(i == nrow(exons)){
                    start_position_long <- marker_codon_end-10
                    end_position_long <- exon_end_position
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(i < nrow(exons)){
                    exon_start_position_after <- exons[i+1,2]
                    exon_end_position_after <- exons[i+1,3]
                    start_position_long <- marker_codon_end-10
                    end_position_long <- exon_end_position
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    number_of_nts <- 10 - (exon_end_position - edit_pos) 
                    exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                    exon_edit <- paste(exon_edit, exon_after, sep = "")
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                }
                else if(left_distance_2 < 10 && right_distance_2 >= 10){
                  if(i == 1){
                    start_position_long <- exon_start_position
                    end_position_long <- edit_pos+10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(i > 1){
                    exon_start_position_before <- exons[i-1,2]
                    exon_end_position_before <- exons[i-1,3]
                    start_position_long <- exon_start_position
                    end_position_long <- edit_pos + 10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    number_of_nts <- 10 - (marker_codon_end - exon_start_position)
                    exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                    exon_edit <- paste( exon_before, exon_edit, sep = "")
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                }
              }
              else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon != "TAA" | marker_codon != "TAG"){
                edit_pos <- exon_start_position  + nt_2
                text <- "ATG, TGG, TAA, TAG and TGA are before the edit and not enough bp after the edit."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(left_distance < 5){
              edit_pos <- exon_start_position  + nt_2
              text <- "Not enough bp before and after the edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
      }
      else if(length_exon > j && counter < 3){
        counter = counter + 1
        j = j+1
      }
      else if(length_exon == j && counter < 3){
        if(counter == 2){
          exon_edit <- gene_seq
          exon_start_position_2 <- exons[i+1,2]
          exon_end_position_2 <- exons[i+1,3]
          exon_seq_2 <- str_sub(gene_seq, exon_start_position_2, exon_end_position_2)
          exon_edit_2_subseq <- exon_seq_2
          chars_from_exon_1 <- str_sub(gene_seq, exon_end_position-1, exon_end_position)
          chars_from_exon_2 <- str_sub(gene_seq, exon_start_position_2, exon_start_position_2)
          chars_combined <- paste(chars_from_exon_1, chars_from_exon_2, sep = "")
          ref <- str_sub(chars_combined, 3, 3)
          var <- get_silent_mutation(chars_combined)
          if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
            edit_pos <- exon_start_position_2
            text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
            t=t+1
            codon_counter = codon_counter + 1
            input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
            input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
            input_list_control[t,3] <- "Synonymous"
            input_list_control[t,4] <- text
          }
          else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
            edit <- var            
            str_sub(exon_edit_2_subseq, 1, 1) <- edit
            str_sub(exon_edit, exon_start_position_2, exon_end_position_2) <- exon_edit_2_subseq
            edit_pos <- exon_start_position_2
            right_distance <- exon_end_position_2 - edit_pos
            if(right_distance >= 3){
              marker_edit_pos <- edit_pos + 3
              marker_codon_start <- edit_pos + 1
              marker_codon_end <- edit_pos + 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                exon_before_start <- exons[i,2]
                exon_before_end <- exons[i,3]
                exon_before_seq <- str_sub(gene_seq, exon_before_end - 10, exon_before_end)
                exon_edit <- str_sub(exon_edit, exon_start_position_2, marker_codon_end +10)
                seq_of_edit <- paste(exon_before_seq, exon_edit, sep = "")
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- seq_of_edit
              }
              else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                edit_pos <- exon_start_position_2
                text <- "Edit in first position of exon from spliced codon shows ATG, TGG, TAA, TAG and TGA after."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(right_distance < 3){
              edit_pos <- exon_start_position_2
              text <- "Following Exon is less then 3 nt so no marker edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
        else if(counter == 1){
          exon_edit <- gene_seq
          exon_start_position_2 <- exons[i+1,2]
          exon_end_position_2 <- exons[i+1,3]
          exon_seq_1 <- str_sub(gene_seq, exon_start_position_2, exon_end_position_2)
          exon_edit_2_subseq <- exon_seq_1
          chars_from_exon_1 <- str_sub(gene_seq, exon_end_position, exon_end_position)
          chars_from_exon_2 <- str_sub(gene_seq, exon_start_position_2, exon_start_position_2+1)
          chars_combined <- paste(chars_from_exon_1, chars_from_exon_2, sep = "")
          ref <- str_sub(chars_combined, 3, 3)
          var <- get_silent_mutation(chars_combined)
          if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
            edit_pos <- exon_start_position_2 + 1
            text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
            t=t+1
            codon_counter = codon_counter + 1
            input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
            input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
            input_list_control[t,3] <- "Synonymous"
            input_list_control[t,4] <- text
          }
          else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
            edit <- var             
            str_sub(exon_edit_2_subseq, 2, 2) <- edit
            str_sub(exon_edit, exon_start_position_2, exon_end_position_2) <- exon_edit_2_subseq
            edit_pos <- exon_start_position_2 + 1
            right_distance <- exon_end_position_2 - edit_pos
            if(right_distance >= 3){
              marker_edit_pos <- edit_pos + 3
              marker_codon_start <- edit_pos + 1
              marker_codon_end <- edit_pos + 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                exon_before_start <- exons[i,2]
                exon_before_end <- exons[i,3]
                exon_before_seq <- str_sub(gene_seq, exon_before_end - 9, exon_before_end)
                exon_edit <- str_sub(exon_edit, exon_start_position_2, marker_codon_end +10)
                seq_of_edit <- paste(exon_before_seq, exon_edit, sep = "")
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- seq_of_edit
              }
              else if(marker_codon == "ATG" |  marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                edit_pos <- exon_start_position_2 + 1
                text <- "Edit in second position of exon from spliced codon shows ATG, TGG, TAA, TAG and TGA after."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(right_distance < 3){
              edit_pos <- exon_start_position_2 + 1
              text <- "Following Exon is less then 3 nt so no marker edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
      }
    }
  }
  else if(exon_seq_2 == exon_seq){
    counter = 1
    codon_end = 1
    for(l in 2:nchar(exon_seq)){
      exon_edit_subseq <- exon_seq
      marker_codon_seq <- gene_seq
      exon_edit <- gene_seq
      if(length_exon >= l && counter == 3){
        counter = 1
        codon_end = codon_end + 3
        nt_1 <- codon_end - 2
        nt_2 <- codon_end - 1
        nt_3 <- codon_end
        codon <- str_sub(exon_seq, nt_1, nt_3)
        ref <- str_sub(codon, 3, 3)
        var <- get_silent_mutation(codon)
        if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
          edit_pos <- exon_start_position  + nt_2
          text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
          t=t+1
          codon_counter = codon_counter + 1
          input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
          input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
          input_list_control[t,3] <- "Synonymous"
          input_list_control[t,4] <- text
        }
        else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
          ref <- str_sub(codon, 3, 3)
          edit <- var 
          str_sub(exon_edit_subseq, nt_3, nt_3) <- edit
          str_sub(exon_edit, exon_start_position, exon_end_position) <- exon_edit_subseq
          edit_pos <- exon_start_position  + nt_2
          right_distance <- exon_end_position - edit_pos
          left_distance <- edit_pos - exon_start_position
          if(right_distance >= 3){
            marker_edit_pos <- edit_pos + 3
            marker_codon_start <- edit_pos + 1
            marker_codon_end <- edit_pos + 3
            marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
            marker_ref <- str_sub(marker_codon, 3, 3)
            if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
              marker_var <- get_silent_mutation(marker_codon)
              marker_edit <- marker_var
              str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
              left_distance_2 <- edit_pos - exon_start_position
              right_distance_2 <- exon_end_position - marker_edit_pos
              if(left_distance_2 >= 10 && right_distance_2 >= 10){
                start_position_long <- edit_pos-10
                end_position_long <- marker_codon_end+10
                exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                codon_counter = codon_counter + 1
                t=t+1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- exon_edit
              }
              else if(left_distance_2 >= 10 && right_distance_2 < 10){
                if(i == nrow(exons)){
                  start_position_long <- edit_pos-10
                  end_position_long <- exon_end_position
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(i < nrow(exons)){
                  exon_start_position_after <- exons[i+1,2]
                  exon_end_position_after <- exons[i+1,3]
                  start_position_long <- edit_pos-10
                  end_position_long <- exon_end_position
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  number_of_nts <- 10 - (exon_end_position - marker_codon_end) 
                  exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                  exon_edit <- paste(exon_edit, exon_after, sep = "")
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
              }
              else if(left_distance_2 < 10 && right_distance_2 >= 10){
                if(i == 1){
                  start_position_long <- exon_start_position
                  end_position_long <- marker_codon_end+10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(i > 1){
                  exon_start_position_before <- exons[i-1,2]
                  exon_end_position_before <- exons[i-1,3]
                  start_position_long <- exon_start_position
                  end_position_long <- marker_codon_end + 10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  number_of_nts <- 10- (edit_pos - exon_start_position)
                  exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                  exon_edit <- paste(exon_before, exon_edit, sep = "")
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
              }
            }
            else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
              if(left_distance >= 5){
                marker_edit_pos <- edit_pos - 3
                marker_codon_start <- edit_pos - 5
                marker_codon_end <- edit_pos - 3
                marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
                marker_ref <- str_sub(marker_codon, 3, 3)
                if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                  marker_var <- get_silent_mutation(marker_codon)
                  marker_edit <- marker_var
                  str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                  left_distance_2 <- marker_edit_pos - exon_start_position
                  right_distance_2 <- exon_end_position - edit_pos
                  if(left_distance_2 >= 10 && right_distance_2 >= 10){
                    start_position_long <- marker_codon_end-10
                    end_position_long <- edit_pos+10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(left_distance_2 >= 10 && right_distance_2 < 10){
                    if(i == nrow(exons)){
                      start_position_long <- marker_codon_end-10
                      end_position_long <- exon_end_position
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                    else if(i < nrow(exons)){
                      exon_start_position_after <- exons[i+1,2]
                      exon_end_position_after <- exons[i+1,3]
                      start_position_long <- marker_codon_end-10
                      end_position_long <- exon_end_position
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      number_of_nts <- 10 - (exon_end_position - edit_pos) 
                      exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                      exon_edit <- paste(exon_edit, exon_after, sep = "")
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                  }
                  else if(left_distance_2 < 10 && right_distance_2 >= 10){
                    if(i == 1){
                      start_position_long <- exon_start_position
                      end_position_long <- edit_pos+10
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                    else if(i > 1){
                      exon_start_position_before <- exons[i-1,2]
                      exon_end_position_before <- exons[i-1,3]
                      start_position_long <- exon_start_position
                      end_position_long <- edit_pos + 10
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      number_of_nts <- 10- (marker_codon_end - exon_start_position)
                      exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                      exon_edit <- paste( exon_before, exon_edit, sep = "")
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                  }
                }
                else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                  edit_pos <- exon_start_position  + nt_2
                  text <- "ATG, TGG, TAA, TAG and TGA are before and after the edit."
                  t=t+1
                  codon_counter = codon_counter + 1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- text
                }
              }
              else if(left_distance < 5){
                edit_pos <- exon_start_position  + nt_2
                text <- "ATG, TGG, TAA, TAG and TGA are after the edit and not enough bp before the edit."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
          }
          else if(right_distance < 3){
            if(left_distance >= 5){
              marker_edit_pos <- edit_pos - 3
              marker_codon_start <- edit_pos - 5
              marker_codon_end <- edit_pos - 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                left_distance_2 <- marker_edit_pos - exon_start_position
                right_distance_2 <- exon_end_position - edit_pos
                if(left_distance_2 >= 10 && right_distance_2 >= 10){
                  start_position_long <- marker_codon_end-10
                  end_position_long <- edit_pos+10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(left_distance_2 >= 10 && right_distance_2 < 10){
                  if(i == nrow(exons)){
                    start_position_long <- marker_codon_end-10
                    end_position_long <- exon_end_position
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(i < nrow(exons)){
                    exon_start_position_after <- exons[i+1,2]
                    exon_end_position_after <- exons[i+1,3]
                    start_position_long <- marker_codon_end-10
                    end_position_long <- exon_end_position
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    number_of_nts <- 10 - (exon_end_position - edit_pos) 
                    exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                    exon_edit <- paste(exon_edit, exon_after, sep = "")
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                }
                else if(left_distance_2 < 10 && right_distance_2 >= 10){
                  if(i == 1){
                    start_position_long <- exon_start_position
                    end_position_long <- edit_pos+10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(i > 1){
                    exon_start_position_before <- exons[i-1,2]
                    exon_end_position_before <- exons[i-1,3]
                    start_position_long <- exon_start_position
                    end_position_long <- edit_pos + 10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    number_of_nts <- 10 - (marker_codon_end - exon_start_position)
                    exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                    exon_edit <- paste( exon_before, exon_edit, sep = "")
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                }
              }
              else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon != "TAA" | marker_codon != "TAG"){
                edit_pos <- exon_start_position  + nt_2
                text <- "ATG, TGG, TAA, TAG and TGA are before the edit and not enough bp after the edit."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(left_distance < 5){
              edit_pos <- exon_start_position  + nt_2
              text <- "Not enough bp before and after the edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
      }
      else if(length_exon > l && counter < 3){
        counter = counter + 1
        l = l+1
      }
      else if(length_exon == l && counter < 3){
        if(counter == 2){
          exon_edit <- gene_seq
          exon_start_position_2 <- exons[i+1,2]
          exon_end_position_2 <- exons[i+1,3]
          exon_seq_2 <- str_sub(gene_seq, exon_start_position_2, exon_end_position_2)
          exon_edit_2_subseq <- exon_seq_2
          chars_from_exon_1 <- str_sub(gene_seq, exon_end_position-1, exon_end_position)
          chars_from_exon_2 <- str_sub(gene_seq, exon_start_position_2, exon_start_position_2)
          chars_combined <- paste(chars_from_exon_1, chars_from_exon_2, sep = "")
          ref <- str_sub(chars_combined, 3, 3)
          var <- get_silent_mutation(chars_combined)
          if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
            edit_pos <- exon_start_position_2
            text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
            t=t+1
            codon_counter = codon_counter + 1
            input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
            input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
            input_list_control[t,3] <- "Synonymous"
            input_list_control[t,4] <- text
          }
          else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
            edit <- var            
            str_sub(exon_edit_2_subseq, 1, 1) <- edit
            str_sub(exon_edit, exon_start_position_2, exon_end_position_2) <- exon_edit_2_subseq
            edit_pos <- exon_start_position_2
            right_distance <- exon_end_position_2 - edit_pos
            if(right_distance >= 3){
              marker_edit_pos <- edit_pos + 3
              marker_codon_start <- edit_pos + 1
              marker_codon_end <- edit_pos + 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                exon_before_start <- exons[i,2]
                exon_before_end <- exons[i,3]
                exon_before_seq <- str_sub(gene_seq, exon_before_end - 10, exon_before_end)
                exon_edit <- str_sub(exon_edit, exon_start_position_2, marker_codon_end +10)
                seq_of_edit <- paste(exon_before_seq, exon_edit, sep = "")
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- seq_of_edit
              }
              else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                edit_pos <- exon_start_position_2
                text <- "Edit in first position of exon from spliced codon shows ATG, TGG, TAA, TAG and TGA after."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(right_distance < 3){
              edit_pos <- exon_start_position_2
              text <- "Following Exon is less then 3 nt so no marker edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
        else if(counter == 1){
          exon_edit <- gene_seq
          exon_start_position_2 <- exons[i+1,2]
          exon_end_position_2 <- exons[i+1,3]
          exon_seq_1 <- str_sub(gene_seq, exon_start_position_2, exon_end_position_2)
          exon_edit_2_subseq <- exon_seq_1
          chars_from_exon_1 <- str_sub(gene_seq, exon_end_position, exon_end_position)
          chars_from_exon_2 <- str_sub(gene_seq, exon_start_position_2, exon_start_position_2+1)
          chars_combined <- paste(chars_from_exon_1, chars_from_exon_2, sep = "")
          ref <- str_sub(chars_combined, 3, 3)
          var <- get_silent_mutation(chars_combined)
          if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
            edit_pos <- exon_start_position_2 + 1
            text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
            t=t+1
            codon_counter = codon_counter + 1
            input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
            input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
            input_list_control[t,3] <- "Synonymous"
            input_list_control[t,4] <- text
          }
          else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
            edit <- var             
            str_sub(exon_edit_2_subseq, 2, 2) <- edit
            str_sub(exon_edit, exon_start_position_2, exon_end_position_2) <- exon_edit_2_subseq
            edit_pos <- exon_start_position_2 + 1
            right_distance <- exon_end_position_2 - edit_pos
            if(right_distance >= 3){
              marker_edit_pos <- edit_pos + 3
              marker_codon_start <- edit_pos + 1
              marker_codon_end <- edit_pos + 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                exon_before_start <- exons[i,2]
                exon_before_end <- exons[i,3]
                exon_before_seq <- str_sub(gene_seq, exon_before_end - 9, exon_before_end)
                exon_edit <- str_sub(exon_edit, exon_start_position_2, marker_codon_end +10)
                seq_of_edit <- paste(exon_before_seq, exon_edit, sep = "")
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- seq_of_edit
              }
              else if(marker_codon == "ATG" |  marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                edit_pos <- exon_start_position_2 + 1
                text <- "Edit in second position of exon from spliced codon shows ATG, TGG, TAA, TAG and TGA after."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(right_distance < 3){
              edit_pos <- exon_start_position_2 + 1
              text <- "Following Exon is less then 3 nt so no marker edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
      }
    }
  }
  else if(exon_seq_1 == exon_seq){
    counter = 1
    codon_end = 2
    for(m in 3:nchar(exon_seq)){
      exon_edit_subseq <- exon_seq
      marker_codon_seq <- gene_seq
      exon_edit <- gene_seq
      if(length_exon >= m && counter == 3){
        counter = 1
        codon_end = codon_end + 3
        nt_1 <- codon_end - 2
        nt_2 <- codon_end - 1
        nt_3 <- codon_end
        codon <- str_sub(exon_seq, nt_1, nt_3)
        ref <- str_sub(codon, 3, 3)
        var <- get_silent_mutation(codon)
        if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
          edit_pos <- exon_start_position  + nt_2
          text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
          t=t+1
          codon_counter = codon_counter + 1
          input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
          input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
          input_list_control[t,3] <- "Synonymous"
          input_list_control[t,4] <- text
        }
        else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
          ref <- str_sub(codon, 3, 3)
          edit <- var 
          str_sub(exon_edit_subseq, nt_3, nt_3) <- edit
          str_sub(exon_edit, exon_start_position, exon_end_position) <- exon_edit_subseq
          edit_pos <- exon_start_position  + nt_2
          right_distance <- exon_end_position - edit_pos
          left_distance <- edit_pos - exon_start_position
          if(right_distance >= 3){
            marker_edit_pos <- edit_pos + 3
            marker_codon_start <- edit_pos + 1
            marker_codon_end <- edit_pos + 3
            marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
            marker_ref <- str_sub(marker_codon, 3, 3)
            if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
              marker_var <- get_silent_mutation(marker_codon)
              marker_edit <- marker_var
              str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
              left_distance_2 <- edit_pos - exon_start_position
              right_distance_2 <- exon_end_position - marker_edit_pos
              if(left_distance_2 >= 10 && right_distance_2 >= 10){
                start_position_long <- edit_pos-10
                end_position_long <- marker_codon_end+10
                exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                codon_counter = codon_counter + 1
                t=t+1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- exon_edit
              }
              else if(left_distance_2 >= 10 && right_distance_2 < 10){
                if(i == nrow(exons)){
                  start_position_long <- edit_pos-10
                  end_position_long <- exon_end_position
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(i < nrow(exons)){
                  exon_start_position_after <- exons[i+1,2]
                  exon_end_position_after <- exons[i+1,3]
                  start_position_long <- edit_pos-10
                  end_position_long <- exon_end_position
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  number_of_nts <- 10 - (exon_end_position - marker_codon_end) 
                  exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                  exon_edit <- paste(exon_edit, exon_after, sep = "")
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
              }
              else if(left_distance_2 < 10 && right_distance_2 >= 10){
                if(i == 1){
                  start_position_long <- exon_start_position
                  end_position_long <- marker_codon_end+10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(i > 1){
                  exon_start_position_before <- exons[i-1,2]
                  exon_end_position_before <- exons[i-1,3]
                  start_position_long <- exon_start_position
                  end_position_long <- marker_codon_end + 10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  number_of_nts <- 10- (edit_pos - exon_start_position)
                  exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                  exon_edit <- paste(exon_before, exon_edit, sep = "")
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
              }
            }
            else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
              if(left_distance >= 5){
                marker_edit_pos <- edit_pos - 3
                marker_codon_start <- edit_pos - 5
                marker_codon_end <- edit_pos - 3
                marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
                marker_ref <- str_sub(marker_codon, 3, 3)
                if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                  marker_var <- get_silent_mutation(marker_codon)
                  marker_edit <- marker_var
                  str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                  left_distance_2 <- marker_edit_pos - exon_start_position
                  right_distance_2 <- exon_end_position - edit_pos
                  if(left_distance_2 >= 10 && right_distance_2 >= 10){
                    start_position_long <- marker_codon_end-10
                    end_position_long <- edit_pos+10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(left_distance_2 >= 10 && right_distance_2 < 10){
                    if(i == nrow(exons)){
                      start_position_long <- marker_codon_end-10
                      end_position_long <- exon_end_position
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                    else if(i < nrow(exons)){
                      exon_start_position_after <- exons[i+1,2]
                      exon_end_position_after <- exons[i+1,3]
                      start_position_long <- marker_codon_end-10
                      end_position_long <- exon_end_position
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      number_of_nts <- 10 - (exon_end_position - edit_pos) 
                      exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                      exon_edit <- paste(exon_edit, exon_after, sep = "")
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                  }
                  else if(left_distance_2 < 10 && right_distance_2 >= 10){
                    if(i == 1){
                      start_position_long <- exon_start_position
                      end_position_long <- edit_pos+10
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                    else if(i > 1){
                      exon_start_position_before <- exons[i-1,2]
                      exon_end_position_before <- exons[i-1,3]
                      start_position_long <- exon_start_position
                      end_position_long <- edit_pos + 10
                      exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                      number_of_nts <- 10- (marker_codon_end - exon_start_position)
                      exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                      exon_edit <- paste( exon_before, exon_edit, sep = "")
                      codon_counter = codon_counter + 1
                      t=t+1
                      input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                      input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                      input_list_control[t,3] <- "Synonymous"
                      input_list_control[t,4] <- exon_edit
                    }
                  }
                }
                else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                  edit_pos <- exon_start_position  + nt_2
                  text <- "ATG, TGG, TAA, TAG and TGA are before and after the edit."
                  t=t+1
                  codon_counter = codon_counter + 1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- text
                }
              }
              else if(left_distance < 5){
                edit_pos <- exon_start_position  + nt_2
                text <- "ATG, TGG, TAA, TAG and TGA are after the edit and not enough bp before the edit."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
          }
          else if(right_distance < 3){
            if(left_distance >= 5){
              marker_edit_pos <- edit_pos - 3
              marker_codon_start <- edit_pos - 5
              marker_codon_end <- edit_pos - 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                left_distance_2 <- marker_edit_pos - exon_start_position
                right_distance_2 <- exon_end_position - edit_pos
                if(left_distance_2 >= 10 && right_distance_2 >= 10){
                  start_position_long <- marker_codon_end-10
                  end_position_long <- edit_pos+10
                  exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                  codon_counter = codon_counter + 1
                  t=t+1
                  input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                  input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                  input_list_control[t,3] <- "Synonymous"
                  input_list_control[t,4] <- exon_edit
                }
                else if(left_distance_2 >= 10 && right_distance_2 < 10){
                  if(i == nrow(exons)){
                    start_position_long <- marker_codon_end-10
                    end_position_long <- exon_end_position
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(i < nrow(exons)){
                    exon_start_position_after <- exons[i+1,2]
                    exon_end_position_after <- exons[i+1,3]
                    start_position_long <- marker_codon_end-10
                    end_position_long <- exon_end_position
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    number_of_nts <- 10 - (exon_end_position - edit_pos) 
                    exon_after <- str_sub(gene_seq, exon_start_position_after, exon_start_position_after + number_of_nts)
                    exon_edit <- paste(exon_edit, exon_after, sep = "")
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                }
                else if(left_distance_2 < 10 && right_distance_2 >= 10){
                  if(i == 1){
                    start_position_long <- exon_start_position
                    end_position_long <- edit_pos+10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                  else if(i > 1){
                    exon_start_position_before <- exons[i-1,2]
                    exon_end_position_before <- exons[i-1,3]
                    start_position_long <- exon_start_position
                    end_position_long <- edit_pos + 10
                    exon_edit <- substr(exon_edit, start_position_long, end_position_long)
                    number_of_nts <- 10 - (marker_codon_end - exon_start_position)
                    exon_before <- str_sub(gene_seq, exon_end_position_before - number_of_nts, exon_end_position_before)
                    exon_edit <- paste( exon_before, exon_edit, sep = "")
                    codon_counter = codon_counter + 1
                    t=t+1
                    input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                    input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                    input_list_control[t,3] <- "Synonymous"
                    input_list_control[t,4] <- exon_edit
                  }
                }
              }
              else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon != "TAA" | marker_codon != "TAG"){
                edit_pos <- exon_start_position  + nt_2
                text <- "ATG, TGG, TAA, TAG and TGA are before the edit and not enough bp after the edit."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(left_distance < 5){
              edit_pos <- exon_start_position  + nt_2
              text <- "Not enough bp before and after the edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
      }
      else if(length_exon > m && counter < 3){
        counter = counter + 1
        m = m+1
      }
      else if(length_exon == m && counter < 3){
        if(counter == 2){
          exon_edit <- gene_seq
          exon_start_position_2 <- exons[i+1,2]
          exon_end_position_2 <- exons[i+1,3]
          exon_seq_2 <- str_sub(gene_seq, exon_start_position_2, exon_end_position_2)
          exon_edit_2_subseq <- exon_seq_2
          chars_from_exon_1 <- str_sub(gene_seq, exon_end_position-1, exon_end_position)
          chars_from_exon_2 <- str_sub(gene_seq, exon_start_position_2, exon_start_position_2)
          chars_combined <- paste(chars_from_exon_1, chars_from_exon_2, sep = "")
          ref <- str_sub(chars_combined, 3, 3)
          var <- get_silent_mutation(chars_combined)
          if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
            edit_pos <- exon_start_position_2
            text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
            t=t+1
            codon_counter = codon_counter + 1
            input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
            input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
            input_list_control[t,3] <- "Synonymous"
            input_list_control[t,4] <- text
          }
          else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
            edit <- var            
            str_sub(exon_edit_2_subseq, 1, 1) <- edit
            str_sub(exon_edit, exon_start_position_2, exon_end_position_2) <- exon_edit_2_subseq
            edit_pos <- exon_start_position_2
            right_distance <- exon_end_position_2 - edit_pos
            if(right_distance >= 3){
              marker_edit_pos <- edit_pos + 3
              marker_codon_start <- edit_pos + 1
              marker_codon_end <- edit_pos + 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                exon_before_start <- exons[i,2]
                exon_before_end <- exons[i,3]
                exon_before_seq <- str_sub(gene_seq, exon_before_end - 10, exon_before_end)
                exon_edit <- str_sub(exon_edit, exon_start_position_2, marker_codon_end +10)
                seq_of_edit <- paste(exon_before_seq, exon_edit, sep = "")
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- seq_of_edit
              }
              else if(marker_codon == "ATG" | marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                edit_pos <- exon_start_position_2
                text <- "Edit in first position of exon from spliced codon shows ATG, TGG, TAA, TAG and TGA after."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(right_distance < 3){
              edit_pos <- exon_start_position_2
              text <- "Following Exon is less then 3 nt so no marker edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
        else if(counter == 1){
          exon_edit <- gene_seq
          exon_start_position_2 <- exons[i+1,2]
          exon_end_position_2 <- exons[i+1,3]
          exon_seq_1 <- str_sub(gene_seq, exon_start_position_2, exon_end_position_2)
          exon_edit_2_subseq <- exon_seq_1
          chars_from_exon_1 <- str_sub(gene_seq, exon_end_position, exon_end_position)
          chars_from_exon_2 <- str_sub(gene_seq, exon_start_position_2, exon_start_position_2+1)
          chars_combined <- paste(chars_from_exon_1, chars_from_exon_2, sep = "")
          ref <- str_sub(chars_combined, 3, 3)
          var <- get_silent_mutation(chars_combined)
          if(var == "ATG" | var == "TGG" | var == "TGA" | var == "TAA" | var == "TAG"){
            edit_pos <- exon_start_position_2 + 1
            text <- "Edit codon is ATG, TGG, TAA, TAG or TGA."
            t=t+1
            codon_counter = codon_counter + 1
            input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
            input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
            input_list_control[t,3] <- "Synonymous"
            input_list_control[t,4] <- text
          }
          else if(var != "ATG" && var != "TGG" && var != "TGA" && var != "TAA" && var != "TAG"){
            edit <- var             
            str_sub(exon_edit_2_subseq, 2, 2) <- edit
            str_sub(exon_edit, exon_start_position_2, exon_end_position_2) <- exon_edit_2_subseq
            edit_pos <- exon_start_position_2 + 1
            right_distance <- exon_end_position_2 - edit_pos
            if(right_distance >= 3){
              marker_edit_pos <- edit_pos + 3
              marker_codon_start <- edit_pos + 1
              marker_codon_end <- edit_pos + 3
              marker_codon <- str_sub(marker_codon_seq, marker_codon_start, marker_codon_end)
              marker_ref <- str_sub(marker_codon, 3, 3)
              if(marker_codon != "ATG" && marker_codon != "TGG" && marker_codon != "TGA" && marker_codon != "TAA" && marker_codon != "TAG"){
                marker_var <- get_silent_mutation(marker_codon)
                marker_edit <- marker_var
                str_sub(exon_edit, marker_edit_pos, marker_edit_pos) <- marker_edit
                exon_before_start <- exons[i,2]
                exon_before_end <- exons[i,3]
                exon_before_seq <- str_sub(gene_seq, exon_before_end - 9, exon_before_end)
                exon_edit <- str_sub(exon_edit, exon_start_position_2, marker_codon_end +10)
                seq_of_edit <- paste(exon_before_seq, exon_edit, sep = "")
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- seq_of_edit
              }
              else if(marker_codon == "ATG" |  marker_codon == "TGG" | marker_codon == "TGA" | marker_codon == "TAA" | marker_codon == "TAG"){
                edit_pos <- exon_start_position_2 + 1
                text <- "Edit in second position of exon from spliced codon shows ATG, TGG, TAA, TAG and TGA after."
                t=t+1
                codon_counter = codon_counter + 1
                input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
                input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
                input_list_control[t,3] <- "Synonymous"
                input_list_control[t,4] <- text
              }
            }
            else if(right_distance < 3){
              edit_pos <- exon_start_position_2 + 1
              text <- "Following Exon is less then 3 nt so no marker edit."
              t=t+1
              codon_counter = codon_counter + 1
              input_list_control[t,1] <- paste0("Synonymous_edit_", edit_pos)
              input_list_control[t,2] <- paste0("codon_pos_", codon_counter)
              input_list_control[t,3] <- "Synonymous"
              input_list_control[t,4] <- text
            }
          }
        }
      }
    }
  }
}

input_list_control_split <- as.data.frame(str_split_fixed(input_list_control$InputID, "_", 3))
input_list_control <- cbind(input_list_control, input_list_control_split$V3)
colnames(input_list_control)[5] <- "edit_pos"
input_list_control$edit_pos <- as.numeric(as.character(input_list_control$edit_pos))


exon27_control_list <- data.frame()
j=1
for (i in 1:nrow(input_list_control)){
  edit_pos <- input_list_control[i,5]
  exon_27_pos_start <- exons[27,2]
  exon_27_pos_end <- exons[27,3]
  if(edit_pos >= exon_27_pos_start && edit_pos <= exon_27_pos_end){
    exon27_control_list[j,1] <- input_list_control[i,1]
    j = j+1
  }
  else if(edit_pos < exon_27_pos_start | edit_pos > exon_27_pos_end){
    i = i+1
  }
} 

colnames(exon27_control_list) <- "InputID"
exon27_control <- merge(input_list_control, exon27_control_list, by = "InputID", all = FALSE)

write_xlsx(exon27_control, "Searching_sequences_synonymous_edits_cDNA.xlsx") 

