#!/bin/bash
#set -o errexit
#set -o pipefail

if [ $# -eq 0 ]; then
echo >&2 "
$(basename $0) - Filtering FASTQ files to include sequences within a certain length range
USAGE: $(basename $0) -i <FASTQ file input> -o <output file> [OPTIONS]
-i    FASTQ file to be filtered                                                       [required]
-o    FASTQ output                                                [required]
NOTES:
"
exit 1
fi

while getopts "i:o:" op
do
    case "$op" in
        i)  input="$OPTARG";;
        o)  output="$OPTARG";;
        \?) exit 1;;
    esac
done

awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 15 && length(seq) <= 70) {print header, seq, qheader, qseq}}' < ${input} > ${output}

