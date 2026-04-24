#!/bin/bash

# gunzip fastq files
for file in *.fastq.gz *.fq.gz; do
    if [[ -e "$file" && ! -e "${file%.gz}" ]]; then
        echo "Unzip $file..."
        gunzip "$file"
    fi
done

# check if we have to concatenate the reads
read1_files=($(ls *_R1_*.fq *_R1_*.fastq *_1.fq *_1.fastq 2>/dev/null))
read2_files=($(ls *_R2_*.fq *_R2_*.fastq *_2.fq *_2.fastq 2>/dev/null))

# function to extract the prefix
extract_prefix() {
  local filename="$1"
  
  # For MiniSeq Data: remove everything after first "_S" (inklusive)
  if [[ "$filename" == *"_S"* ]]; then
    echo "$filename" | sed -E 's/_S[0-9]+_L[0-9]+_R?[12]_[0-9]+\.?(fastq|fq)?$//'
  
  # For NovoSeq data: remove specific pattern
  elif [[ "$filename" == *"_"* ]]; then
    # Remove pattern EKDL240031543-1A_HK3FNDSXC_L4 or similar
    echo "$filename" | sed -E 's/_[A-Z0-9-]+_[A-Z0-9]+_L[0-9]+_[12]\.(fastq|fq)$//'
  else
    echo "$filename"
  fi
}

# concatenate reads for read 1 and read 2
concat_files() {
  local prefix="$1"
  local read1_files=($(ls *_R1_*.fq *_R1_*.fastq *_1.fq *_1.fastq 2>/dev/null))
  local read2_files=($(ls *_R2_*.fq *_R2_*.fastq *_2.fq *_2.fastq 2>/dev/null))

  echo "Read 1 data found: ${read1_files[@]}"
  echo "Read 2 data found: ${read2_files[@]}"

  if [ ${#read1_files[@]} -gt 1 ]; then
    cat "${read1_files[@]}" > "${prefix}_R1_concatenate.fastq"
    echo "Read 1 data are concatenated: ${prefix}_R1_concatenate.fastq"
  elif [ ${#read1_files[@]} -eq 1 ]; then
    echo "Only one read 1 exists."
  else
    echo "Zero read 1 with $prefix could be found."
    return 1
  fi

  if [ ${#read2_files[@]} -gt 1 ]; then
    cat "${read2_files[@]}" > "${prefix}_R2_concatenate.fastq"
    echo "Read 2 data are concatenated: ${prefix}_R2_concatenate.fastq"
  elif [ ${#read2_files[@]} -eq 1 ]; then
    echo "Only one read 2 exists."
  else
    echo "Zero read 2 with $prefix could be found."
    return 1
  fi

  return 0
}

# Main workflow
main() {
  if [[ ${#read1_files[@]} -eq 0 ]]; then
    echo "No read 1 data could be found."
    exit 1
  fi

  prefix=$(extract_prefix "${read1_files[0]}")
  concat_files "$prefix"

  # Check if concatenated data exists
  if [[ -e "${prefix}_R1_concatenate.fastq" && -e "${prefix}_R2_concatenate.fastq" ]]; then
    echo "Run FLASH: ${prefix}_R1_concatenate.fastq and ${prefix}_R2_concatenate.fastq"
    flash "${prefix}_R1_concatenate.fastq" "${prefix}_R2_concatenate.fastq" -o "$prefix"
  else
    # Names for NovaSeq oder MiniSeq
    novaseq_r1=($(ls "${prefix}"*_1.fq "${prefix}"*_1.fastq 2>/dev/null))
    miniseq_r1=($(ls "${prefix}"*_R1_[0-9]*.fastq "${prefix}"*_R1_[0-9]*.fq 2>/dev/null))
    novaseq_r2=($(ls "${prefix}"*_2.fq "${prefix}"*_2.fastq 2>/dev/null))
    miniseq_r2=($(ls "${prefix}"*_R2_[0-9]*.fastq "${prefix}"*_R2_[0-9]*.fq 2>/dev/null))

    if [[ -n "${novaseq_r1[0]}" && -n "${novaseq_r2[0]}" ]]; then
      echo "Run FLASH with NovaSeq data: ${novaseq_r1[0]} and ${novaseq_r2[0]}"
      flash "${novaseq_r1[0]}" "${novaseq_r2[0]}" -o "$prefix"
    elif [[ -n "${miniseq_r1[0]}" && -n "${miniseq_r2[0]}" ]]; then
      echo "Run FLASH with MiniSeq data: ${miniseq_r1[0]} and ${miniseq_r2[0]}"
      flash "${miniseq_r1[0]}" "${miniseq_r2[0]}" -o "$prefix"
    else
      echo "Could not run FLASH because of missing data"
    fi
  fi
}

main

# gzip .fq and .fastq files
echo "gzip all .fq und .fastq files..."
for file in *.fq *.fastq; do
  if [ -f "$file" ]; then
    gzip "$file"
  fi
done


# start python script by masking low quality base and search gRNAs
echo "Start python script."
python3 mask_low_quality_bases.py "${prefix}.extendedFrags.fastq.gz" "${prefix}_masked.fastq.gz"


# change masked fastq into txt
gunzip *_masked.fastq.gz
for file in *_masked.fastq; do 
  sed -n '/^@/{n;p;}' "$file" > ${file/%.fastq/.txt}; 
done

awk '{ print $2 }' Search_sequences_gDNA.txt > sequences.txt

for i in *_masked.txt; do
  echo "$i"
  dir=${i%%_masked.txt}
  dir=${dir##*/}
  mkdir "${dir%.txt}"
  wc -l $i >> "${dir%.txt}"/results_${dir}.txt

  while read seq; do
    echo "$seq"
    ggrep -i "$seq" $i > "${dir%.txt}"/Edit_${seq}.txt
    cd ${dir%.txt}
    wc -l Edit_${seq}.txt >> results_${dir}.txt
    rm Edit_${seq}.txt
    cd ..
  done < sequences.txt
done

rm sequences.txt

# gzip .fq and .fastq files
echo "gzip all .fq und .fastq files..."
for file in *.fq *.fastq; do
  if [ -f "$file" ]; then
    gzip "$file"
  fi
done
