rm(list=ls())
library(readxl)
library(writexl)
library(stringi)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)


setwd("~/Documents/Projekte/Prime_Editing/PE_Paper_20241126/Zenodo/Design_pegRNAs_cBioPortal_mutations")

# Input sequence of the gene and the lists from cBioPortal 
gene_seq <- read_file("Chr17_full_sequence.txt")
nchar(gene_seq)

mutations <- read_tsv("cBioPortal_27Jan2022.tsv" )
#remove all rows that contain "fusion" and "Splice-Site" mutation types
mutations <- as.data.frame(filter(mutations, mutations$`Mutation Type`!="fusion"))
mutations <- as.data.frame(filter(mutations, mutations$`Mutation Type`!="Splice_Site"))
mutations <- as.data.frame(filter(mutations, mutations$`Mutation Type`!="Splice_Region"))
mutations <- as.data.frame(filter(mutations, !is.na(mutations$`Variant Type`)))

# combine names of Sample ID with Variant Type to generate unique Sample ID's and remove duplicates
mutations_unique <- mutations[!duplicated(mutations$`Protein Change`), ]
mutations_unique$`Sample ID` <- paste(mutations_unique$`Sample ID`, mutations_unique$`Protein Change`, sep = "_")
mutations_unique <- mutations_unique[!duplicated(mutations_unique$`Sample ID`), ]
nr_of_mutations <- table(mutations_unique$`Mutation Type`)
nr_of_mutations <- as.data.frame(nr_of_mutations)
colnames(nr_of_mutations) <- c("Mutation_Type", "count")
#write_xlsx(nr_of_mutations, "Mutation_count.xlsx")


# read in coordinates of exons 
exons <- as.data.frame(read.table("NF1-Exons-positions.txt"))
exons <- exons[-c(1,2),]
colnames(exons) <- c("Exon", "Start", "End")
exons <- as.data.frame(sapply(exons,as.numeric))

## read in exon splice information
exon_splice_information <- read.table("Exon_splice_information.txt", header = TRUE)

# get exon
get_exon_info <- function(which_exon){
  if(is.na(which_exon)){
    for (z in 1:nrow(exons)){
      if (start_position >= exons[z,2] && start_position <= exons[z,3]){
        which_exon <- exons[z,1]
        return(which_exon)
      }
      else{
        z = z+1
      }
    }
  }
  else{
    which_exon <- as.data.frame(str_split(which_exon, "/"))
    which_exon <- which_exon[1,1]
    return(which_exon)
  }
}

# calculate codon position of edit
get_codon_pos_info <- function(exon_splice_info_begin){
  if(exon_splice_info_begin == 3){
    codon_pos <- format((edit_pos - exon_start_pos)/3, nsmall = 2)
    codon_pos <- as.data.frame(unlist(strsplit(codon_pos, "[.]")))
    codon_pos <- codon_pos[2,1]
    codon_pos <- substr(codon_pos, start = 1, stop = 2)
    codon_pos <- as.numeric(codon_pos)
    return(codon_pos)
  }
  else if(exon_splice_info_begin == 2){
    codon_pos <- format((edit_pos - exon_start_pos + 1)/3, nsmall = 2)
    codon_pos <- as.data.frame(unlist(strsplit(codon_pos, "[.]")))
    codon_pos <- codon_pos[2,1]
    codon_pos <- substr(codon_pos, start = 1, stop = 2)
    codon_pos <- as.numeric(codon_pos)
    return(codon_pos)
  }
  else if(exon_splice_info_begin == 1){
    codon_pos <- format((edit_pos - exon_start_pos - 1)/3, nsmall = 2)
    codon_pos <- as.data.frame(unlist(strsplit(codon_pos, "[.]")))
    codon_pos <- codon_pos[2,1]
    codon_pos <- substr(codon_pos, start = 1, stop = 2)
    codon_pos <- as.numeric(codon_pos)
    return(codon_pos)
  }
}

# calculate cododn position of edit start
get_codon_pos_info_start <- function(exon_splice_info_begin){
  if(exon_splice_info_begin == 3){
    codon_pos <- format((start_position - exon_start_pos)/3, nsmall = 2)
    codon_pos <- as.data.frame(unlist(strsplit(codon_pos, "[.]")))
    codon_pos <- codon_pos[2,1]
    codon_pos <- substr(codon_pos, start = 1, stop = 2)
    codon_pos <- as.numeric(codon_pos)
    return(codon_pos)
  }
  else if(exon_splice_info_begin == 2){
    codon_pos <- format((start_position - exon_start_pos + 1)/3, nsmall = 2)
    codon_pos <- as.data.frame(unlist(strsplit(codon_pos, "[.]")))
    codon_pos <- codon_pos[2,1]
    codon_pos <- substr(codon_pos, start = 1, stop = 2)
    codon_pos <- as.numeric(codon_pos)
    return(codon_pos)
  }
  else if(exon_splice_info_begin == 1){
    codon_pos <- format((start_position - exon_start_pos - 1)/3, nsmall = 2)
    codon_pos <- as.data.frame(unlist(strsplit(codon_pos, "[.]")))
    codon_pos <- codon_pos[2,1]
    codon_pos <- substr(codon_pos, start = 1, stop = 2)
    codon_pos <- as.numeric(codon_pos)
    return(codon_pos)
  }
}

##################################################################################################################

## create input list for WT sequences based on cBioPortal mutations
input_list <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(input_list) <- c("InputID", "MutationType", "InputSequence")   
gene_seq_edit <- gene_seq 
for (i in 1:nrow(mutations_unique)){
  gene_seq_edit <- gene_seq 
  exon_seq_marker <- gene_seq
  if(mutations_unique[i,9] == "SNP"){
    inputID <- mutations_unique[i,2]
    input_list[i,1] <- inputID
    mutation_type <- mutations_unique[i,8]
    start_position <- mutations_unique[i,16]
    end_position <- mutations_unique[i,17]
    which_exon <- mutations_unique[i,29]
    which_exon <- get_exon_info(which_exon)
    # find out on which codon position edit is and include marker edit
    exon <- exons[which_exon,]
    exon_start_pos <- exon[,2]
    exon_end_pos <- exon[,3]
    exon_splice_info_begin <- exon_splice_information[which_exon,2]
    exon_splice_info_end <- exon_splice_information[which_exon,3]
    edit_pos <- end_position
    right_distance <- exon_end_pos - edit_pos
    left_distance <- edit_pos - exon_start_pos
    codon_pos <- get_codon_pos_info(exon_splice_info_begin)
    if(right_distance >= 4 && codon_pos == 33){
      start_pos_marker_edit <- end_position + 4
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 2, end_position + 4)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(left_distance >= 4){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, end_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 4 && codon_pos == 33){
      if(left_distance >= 4){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, end_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "No full codon on right side from edit and only codon with ATG, TGG, TAA, TAG or TGA before the edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(left_distance < 4){
        text <- "No full codon on the left and right site of edit."
        input_list[i,2] <- mutation_type
        input_list[i,3] <- text
      }
    }
    else if(right_distance >= 3 && codon_pos == 66){
      start_pos_marker_edit <- end_position + 3
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 1, end_position + 3)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, end_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if (right_distance < 3 && codon_pos == 66){
      if(left_distance >= 5){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, end_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "No full codon on right side from edit and only codon with ATG, TGG, TAA, TAG or TGA before the edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(left_distance < 5){
        text <- "No full codon on the left and right site of edit."
        input_list[i,2] <- mutation_type
        input_list[i,3] <- text
      }
    }
    else if(right_distance >= 5 && codon_pos == 0){
      start_pos_marker_edit <- end_position + 5
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 3, end_position + 5)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(left_distance >= 3){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, end_position - 1)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 5 && codon_pos == 0){
      if(left_distance >= 3){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, end_position - 1)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "No full codon on right side from edit and only codon with ATG, TGG, TAA, TAG or TGA before the edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(left_distance < 3){
        text <- "No full codon on the left and right site of edit."
        input_list[i,2] <- mutation_type
        input_list[i,3] <- text
      }
    }
  }
  else if(mutations_unique[i,9] == "DEL"){
    inputID <- mutations_unique[i,2]
    input_list[i,1] <- inputID
    mutation_type <- mutations_unique[i,8]
    start_position <- mutations_unique[i,16]
    end_position <- mutations_unique[i,17]
    start_position_long <- start_position-10
    end_position_long <- end_position+10
    gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
    chars <- "-"
    if(grepl(chars, gene_seq_edit, fixed = TRUE)){
      gene_seq_edit <- str_remove_all(gene_seq_edit, "-")
      input_list[i,2] <- gene_seq_edit
    }
    else{
      input_list[i,2] <- gene_seq_edit
    }
    input_list[i,2] <- mutation_type
    input_list[i,3] <- gene_seq_edit
  }
  else if(mutations_unique[i,9] == "INS"){
    inputID <- mutations_unique[i,2]
    input_list[i,1] <- inputID
    mutation_type <- mutations_unique[i,8]
    start_position <- mutations_unique[i,16]
    end_position <- mutations_unique[i,17]
    start_position_long <- start_position-10
    end_position_long <- end_position+10
    gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
    input_list[i,2] <- mutation_type
    input_list[i,3] <- gene_seq_edit
  }
  else if(mutations_unique[i,9] == "DNP"){
    inputID <- mutations_unique[i,2]
    input_list[i,1] <- inputID
    start_position <- mutations_unique[i,16]
    end_position <- mutations_unique[i,17]
    which_exon <- mutations_unique[i,29]
    which_exon <- get_exon_info(which_exon)
    # find out on which codon position edit is and include marker edit
    exon <- exons[which_exon,]
    exon_start_pos <- exon[,2]
    exon_end_pos <- exon[,3]
    exon_splice_info_begin <- exon_splice_information[which_exon,2]
    exon_splice_info_end <- exon_splice_information[which_exon,3]
    edit_pos <- end_position
    right_distance <- exon_end_pos - edit_pos
    left_distance <- edit_pos - exon_start_pos
    codon_pos <- get_codon_pos_info(exon_splice_info_begin)
    if(right_distance >= 4 && codon_pos == 33){
      start_pos_marker_edit <- end_position + 4
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 2, end_position + 4)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(left_distance >= 4){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, end_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 4 && codon_pos == 33){
      if(left_distance >= 4){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, end_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "No full codon on right side from edit and only codon with ATG, TGG, TAA, TAG or TGA before the edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(left_distance < 4){
        text <- "No full codon on the left and right site of edit."
        input_list[i,2] <- mutation_type
        input_list[i,3] <- text
      }
    }
    else if(right_distance >= 3 && codon_pos == 66){
      start_pos_marker_edit <- end_position + 3
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 1, end_position + 3)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(left_distance >= 5){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, end_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
            
          }
        }
      }
    }
    else if (right_distance < 3 && codon_pos == 66){
      if(left_distance >= 5){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, end_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "No full codon on right side from edit and only codon with ATG, TGG, TAA, TAG or TGA before the edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(left_distance < 5){
        text <- "No full codon on the left and right site of edit."
        input_list[i,2] <- mutation_type
        input_list[i,3] <- text
      }
    }
    else if(right_distance >= 5 && codon_pos == 0){
      start_pos_marker_edit <- end_position + 5
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 3, end_position + 5)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, end_position - 4)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 5 && codon_pos == 0){
      if(left_distance >= 4){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, end_position - 4)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "No full codon on right side from edit and only codon with ATG, TGG, TAA, TAG or TGA before the edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(left_distance < 4){
        text <- "No full codon on the left and right site of edit."
        input_list[i,2] <- mutation_type
        input_list[i,3] <- text
      }
    }
  }
  else if(mutations_unique[i,9] == "ONP"){
    inputID <- mutations_unique[i,2]
    input_list[i,1] <- inputID
    mutation_type <- mutations_unique[i,8]
    start_position <- mutations_unique[i,16]
    end_position <- mutations_unique[i,17]
    which_exon <- mutations_unique[i,29]
    which_exon <- get_exon_info(which_exon)
    # find out on which codon position edit is and include marker edit
    exon <- exons[which_exon,]
    exon_start_pos <- exon[,2]
    exon_end_pos <- exon[,3]
    exon_splice_info_begin <- exon_splice_information[which_exon,2]
    exon_splice_info_end <- exon_splice_information[which_exon,3]
    edit_pos <- end_position
    right_distance <- exon_end_pos - edit_pos
    left_distance <- edit_pos - exon_start_pos
    codon_pos <- get_codon_pos_info(exon_splice_info_begin)
    codon_pos_start <- get_codon_pos_info_start(exon_splice_info_begin)
    if(right_distance >= 4 && codon_pos == 33){
      start_pos_marker_edit <- end_position + 4
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 2, end_position + 4)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(codon_pos_start == 0 && left_distance >= 3){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 33 && left_distance >= 4){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 66 && left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 4 && codon_pos == 33){
      if(codon_pos_start == 0 && left_distance >= 3){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 33 && left_distance >= 4){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 66 && left_distance >= 5){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
    }
    else if(right_distance >= 3 && codon_pos == 66){
      start_pos_marker_edit <- end_position + 3
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 1, end_position + 3)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(codon_pos_start == 0 && left_distance >= 3){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 33 && left_distance >= 4){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 66 && left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 3 && codon_pos == 66){
      if(codon_pos_start == 0 && left_distance >= 3){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 33 && left_distance >= 4){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 66 && left_distance >= 5){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
    }
    else if(right_distance >= 5 && codon_pos == 0){
      start_pos_marker_edit <- end_position + 5
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 3, end_position + 5)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(codon_pos_start == 0 && left_distance >= 3){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 33 && left_distance >= 4){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 66 && left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 5 && codon_pos == 0){
      if(codon_pos_start == 0 && left_distance >= 3){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 33 && left_distance >= 4){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 66 && left_distance >= 5){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
    }
  }
  else if(mutations_unique[i,9] == "TNP"){
    inputID <- mutations_unique[i,2]
    input_list[i,1] <- inputID
    start_position <- mutations_unique[i,16]
    end_position <- mutations_unique[i,17]
    which_exon <- mutations_unique[i,29]
    which_exon <- get_exon_info(which_exon)
    # find out on which codon position edit is and include marker edit
    exon <- exons[which_exon,]
    exon_start_pos <- exon[,2]
    exon_end_pos <- exon[,3]
    exon_splice_info_begin <- exon_splice_information[which_exon,2]
    exon_splice_info_end <- exon_splice_information[which_exon,3]
    edit_pos <- end_position
    right_distance <- exon_end_pos - edit_pos
    left_distance <- edit_pos - exon_start_pos
    codon_pos <- get_codon_pos_info(exon_splice_info_begin)
    codon_pos_start <- get_codon_pos_info_start(exon_splice_info_begin)
    if(right_distance >= 4 && codon_pos == 33){
      start_pos_marker_edit <- end_position + 4
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 2, end_position + 4)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(codon_pos_start == 0 && left_distance >= 3){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 33 && left_distance >= 4){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            iinput_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 66 && left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 4 && codon_pos == 33){
      if(codon_pos_start == 0 && left_distance >= 3){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 33 && left_distance >= 4){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 66 && left_distance >= 5){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
    }
    else if(right_distance >= 3 && codon_pos == 66){
      start_pos_marker_edit <- end_position + 3
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 1, end_position + 3)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(codon_pos_start == 0 && left_distance >= 3){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 33 && left_distance >= 4){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 66 && left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if (right_distance < 3 && codon_pos == 66){
      if(codon_pos_start == 0 && left_distance >= 3){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 33 && left_distance >= 4){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 66 && left_distance >= 5){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
    }
    else if(right_distance >= 5 && codon_pos == 0){
      start_pos_marker_edit <- end_position + 5
      end_pos_marker_edit <- start_pos_marker_edit
      marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
      marker_full_codon <- substr(exon_seq_marker, end_position + 3, end_position + 5)
      if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
        start_position_long <- start_position-10
        end_position_long <- end_pos_marker_edit+10
        gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
        input_list[i,2] <- mutation_type
        input_list[i,3] <- gene_seq_edit
      }
      else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
        if(codon_pos_start == 0 && left_distance >= 3){
          start_pos_marker_edit <- start_position - 1
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 33 && left_distance >= 4){
          start_pos_marker_edit <- start_position - 2
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
        else if(codon_pos_start == 66 && left_distance >= 5){
          start_pos_marker_edit <- start_position - 3
          end_pos_marker_edit <- start_pos_marker_edit
          marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
          marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
          if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
            start_position_long <- start_pos_marker_edit-10
            end_position_long <- end_position+10
            gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
            input_list[i,2] <- mutation_type
            input_list[i,3] <- gene_seq_edit
          }
          else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
            text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
            input_list[i,2] <- mutation_type
            input_list[i,3] <- text
          }
        }
      }
    }
    else if(right_distance < 5 && codon_pos == 0){
      if(codon_pos_start == 0 && left_distance >= 3){
        start_pos_marker_edit <- start_position - 1
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 3, start_position - 1)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 33 && left_distance >= 4){
        start_pos_marker_edit <- start_position - 2
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 4, start_position - 2)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
      else if(codon_pos_start == 66 && left_distance >= 5){
        start_pos_marker_edit <- start_position - 3
        end_pos_marker_edit <- start_pos_marker_edit
        marker_ref <- substr(exon_seq_marker, start_pos_marker_edit, end_pos_marker_edit)
        marker_full_codon <- substr(exon_seq_marker, start_position - 5, start_position - 3)
        if(marker_full_codon != "ATG" && marker_full_codon != "TGG" && marker_full_codon != "TGA" && marker_full_codon != "TAA" && marker_full_codon != "TAG"){
          start_position_long <- start_pos_marker_edit-10
          end_position_long <- end_position+10
          gene_seq_edit <- substr(gene_seq_edit, start_position_long, end_position_long)
          input_list[i,2] <- mutation_type
          input_list[i,3] <- gene_seq_edit
        }
        else if(marker_full_codon == "ATG" | marker_full_codon == "TGG" | marker_full_codon == "TGA" | marker_full_codon == "TAA" | marker_full_codon == "TAG"){
          text <- "Codon is ATG, TGG, TAA, TAG or TGA before and after edit."
          input_list[i,2] <- mutation_type
          input_list[i,3] <- text
        }
      }
    }
  }
  else{
    i = i+1
  }
}


## get only edits for exon 27
exon27_mutation_list <- data.frame()
j=1
for (i in 1:nrow(mutations_unique)){
  if(mutations_unique[i,16] >= exons[27,2] && mutations_unique[i,17] <= exons[27,3]){
    exon27_mutation_list[j,1] <- mutations_unique[i,2]
    exon27_mutation_list[j,2] <- mutations_unique[i,16]
    j = j+1
  }
  else{
    i = i+1
  }
} 

colnames(exon27_mutation_list) <- c("InputID","edit_pos")
exon27_mutations <- merge(input_list, exon27_mutation_list, by = "InputID")

write_xlsx(exon27_mutations, "Searching_sequences_WT_mutations_gDNA.xlsx") 

