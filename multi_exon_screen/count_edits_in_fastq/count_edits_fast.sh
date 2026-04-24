#!/bin/bash
##set -euo pipefail
shopt -s extglob nullglob

helpFunction()
{
    echo ""
    echo "Usage: $0 -p path -s path -r"
    echo "Run cDNA count on a set of fastq files. This script requires python3 and FLASH (Fast Length Adjustment of SHort reads) to be present in your PATH"
    echo -e "\t-p path to folder containing FASTQ files"
    echo -e "\t-s path to searching sequences.txt file"
    echo -e "\t-r optional: provide if reverse complemented search sequences should also be searched. Default: false"
    echo -e "\t-m optional: minimal overlap required for merging r1 and r0. Default: 10"
    exit 1 # Exit script after printing help
}

do_rev_comp=false
min_overlap=10
while getopts "p:s:m:ra" opt
do
    case "$opt" in
        p ) fastqdir="$OPTARG" ;;
        s ) ssfile="$OPTARG" ;;
        m ) min_overlap="$OPTARG" ;;
        r ) do_rev_comp=true ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "$fastqdir" ] || [ -z "$ssfile" ]
then
    echo "ERROR: Some or all of the parameters are empty";
    helpFunction
fi

# Begin script in case all parameters are correct
echo "Counting edits in $fastqdir using searching sequences from $ssfile..."

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd "$fastqdir"

# check if we have to concatenate the reads
#read1_files=($(ls *_1.fq.gz *_1.fastq.gz 2>/dev/null))
#read2_files=($(ls *_2.fq.gz *_2.fastq.gz 2>/dev/null))
read1_files=( !(*notCombined*)_R1_[0-9]*.fq.gz !(*notCombined*)_R1_[0-9]*.fastq.gz !(*notCombined*)_1.fq.gz !(*notCombined*)_1.fastq.gz ) # 
read2_files=( !(*notCombined*)_R2_[0-9]*.fq.gz !(*notCombined*)_R2_[0-9]*.fastq.gz !(*notCombined*)_2.fq.gz !(*notCombined*)_2.fastq.gz ) # 


# function to extract the prefix
extract_prefix() {
    local filename="$1"
  
    # For MiniSeq Data: remove everything after first "_S" (inklusive)
    if [[ "$filename" == *"_S"* ]]; then
        echo "$filename" | sed -E 's/_S[0-9]+_L[0-9]+_R?[12]_[0-9]+\.?(fastq\.gz|fq\.gz)?$//'
  
    # For NovoSeq data: remove specific pattern
    elif [[ "$filename" == *"_"* ]]; then
        # Remove pattern EKDL240031543-1A_HK3FNDSXC_L4 or similar
        echo "$filename" | sed -E 's/_[A-Z0-9-]+_[A-Z0-9]+_L[0-9]+_[12]\.(fastq\.gz|fq\.gz)$//'
    else
        echo "$filename"
    fi
}

# concatenate reads for read 1 and read 2
concat_files() {
    local prefix="$1"
    local -n read1_files_ref="$2"
    local -n read2_files_ref="$3"

    #read1_files=($(printf '%s\n' "${read1_files[@]}" | sed '/notCombined/d'))
    #read2_files=($(printf '%s\n' "${read2_files[@]}" | sed '/notCombined/d'))

    echo "Read 1 data found: ${read1_files_ref[*]}"
    echo "Read 2 data found: ${read2_files_ref[*]}"

    if [ ${#read1_files_ref[@]} -gt 1 ]; then
        echo "Concatenating read 1 data to file: ${prefix}_R1_concatenate.fastq.gz"
        cat "${read1_files_ref[@]}" > "${prefix}_R1_concatenate.fastq.gz"
    elif [ ${#read1_files_ref[@]} -eq 1 ]; then
        echo "one read 1 file found"
    else
        echo "ERROR: no read 1 files found"
        return 1
    fi

    if [ ${#read2_files_ref[@]} -gt 1 ]; then
        echo "Concatenating read 2 data to file: ${prefix}_R2_concatenate.fastq.gz"
        cat "${read2_files_ref[@]}" > "${prefix}_R2_concatenate.fastq.gz"
    elif [ ${#read2_files_ref[@]} -eq 1 ]; then
        echo "one read 2 file found"
    else
        echo "ERROR: no read 2 files found"
        return 1
    fi

    return 0
}

# start script
if [[ ${#read1_files[@]} -eq 0 ]]; then
    echo "ERROR: no read 1 data found."
    exit 1
fi


prefix=$(extract_prefix "${read1_files[0]}")
echo "extracted prefix: $prefix"


# concatename multiple r1 / r2 files in one
if [[ -e "${prefix}_R1_concatenate.fastq.gz" && -e "${prefix}_R2_concatenate.fastq.gz" ]]; then
    echo "SKIP CONCATENATION: found concatenated files"
else
    concat_files "$prefix" read1_files read2_files
fi

# check if concatenated data exists - if yes combine it with flash if not there is only one r1 and r2 fastq. so use it directly
read_1=""
read_2=""
if [[ -e "${prefix}_R1_concatenate.fastq.gz" && -e "${prefix}_R2_concatenate.fastq.gz" ]]; then
    read_1="${prefix}_R1_concatenate.fastq.gz"
    read_2="${prefix}_R2_concatenate.fastq.gz"
elif [[ ${#read1_files[@]} -eq 1 && ${#read2_files[@]} -eq 1 ]]; then
    read_1="${read1_files[0]}"
    read_2="${read2_files[0]}"
fi

# check if flash output exists already
if [[ -e "${prefix}.extendedFrags.fastq.gz" && -e "${prefix}.notCombined_1.fastq.gz" && -e "${prefix}.notCombined_2.fastq.gz" ]]; then
    echo "SKIP MERGE: found flash output"
else
    # run flash
    if [[ -n "$read_1" && -n "$read_2" ]]; then
        echo "running FLASH on $read_1 and $read_2"
        flash "$read_1" "$read_2" -m $min_overlap -z -q -o "$prefix"
    else
        echo "ERROR: could not run FLASH because of missing or ambiguous data"
        exit 1
    fi
fi


# merge output of FLASH
concatfn="${prefix}.extendedFrags.fastq.gz"
#concatfn="${prefix}_concatenated_flash.fastq.gz"
#if [[ -e $concatfn ]]; then
#    echo "SKIP FLASH CONCATENATION: concatenated flash output found"
#else
#    echo "concatenating FLASH output"
#    seqtk seq -r "${prefix}.notCombined_2.fastq.gz" | gzip > "${prefix}.notCombined_2_rv.fastq.gz" # reverse complement read 2
#    cat "${prefix}.extendedFrags.fastq.gz" "${prefix}.notCombined_1.fastq.gz" "${prefix}.notCombined_2_rv.fastq.gz" > $concatfn
#fi


# mask low quality bases
concatfn_masked="${prefix}_concatenated_flash_masked.fastq.gz"
if [[ -e $concatfn_masked ]]; then
    echo "SKIP MASKING: file already exists"
else
    echo "masking low quality reads"
    python3 $SCRIPT_DIR/mask_low_quality_bases.py "$concatfn" "$concatfn_masked"
    #python3 $SCRIPT_DIR/process_counts_nomask.py "$concatfn" $ssfile "${prefix}_counts.tsv"
fi


# count forward search sequences
counting_results_foward="${prefix}_results.txt"
if [[ -e $counting_results_foward ]]; then
    echo "SKIP COUNTING: forward count file exists"
else
    # start counting
    echo "counting forward strand"
    python3 $SCRIPT_DIR/count_ss.py "$ssfile" "$concatfn_masked" "$counting_results_foward"
fi


# count reverse sequences
counting_results_reverse="NONE"
if [[ "$do_rev_comp" == "true" ]]; then
    # write reverse complement search sequences if they do not exist
    ssfile_rv=$ssfile.rv
    if [[ -e $ssfile_rv ]]; then
        echo "SKIP WRITING REVERSE COMPLEMENT SEARCH SEQUENCES: file exists"
    else
        python3 $SCRIPT_DIR/merge_output.py "$ssfile" "$ssfile_rv"
    fi

    # count reverse search sequences
    counting_results_reverse="${prefix}_results_RC.txt"
    if [[ -e $counting_results_reverse ]]; then
        echo "SKIP COUNTING: reverse count file exists"
    else
        # start counting
        echo "counting reverse strand"
        python3 $SCRIPT_DIR/count_ss.py "$ssfile_rv" "$concatfn_masked" "$counting_results_reverse"
    fi
fi


counting_results="${prefix}_results_total.tsv"
if [[ -e $counting_results ]]; then
    echo "SKIP MERGE: merged output exists"
else
    # merge count output - if we have forward and reverse counts then counting_results_reverse is not "NONE"
    echo "mergeing counting result"
    python3 $SCRIPT_DIR/merge_output.py "$counting_results_foward" "$counting_results_reverse" "$ssfile" "$counting_results"

    # calculate the total number of reads
    #"${prefix}.extendedFrags.fastq.gz" "${prefix}.notCombined_1.fastq.gz" "${prefix}.notCombined_2.fastq.gz"
    merged=$(($(zcat "${prefix}.extendedFrags.fastq.gz" | wc -l) / 4))
    unmerged1=$(($(zcat "${prefix}.notCombined_1.fastq.gz" | wc -l) / 4))
    unmerged2=$(($(zcat "${prefix}.notCombined_2.fastq.gz" | wc -l) / 4))

    sed -i "1s/^/#total_number_of_reads_notCombined_2: $unmerged2\n/" "$counting_results"
    sed -i "1s/^/#total_number_of_reads_notCombined_1: $unmerged1\n/" "$counting_results"
    sed -i "1s/^/#total_number_of_reads_combined: $merged\n/" "$counting_results"
fi


echo "~~ DONE ~~"

