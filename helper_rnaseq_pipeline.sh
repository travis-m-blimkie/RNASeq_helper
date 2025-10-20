#!/bin/bash
set -e
p_flag=0
VERSION="0.1.0"

# Preprocess long options into short ones
for arg in "$@"; do
    shift
    case "$arg" in
        --file)       set -- "$@" -f ;;
        --directory)  set -- "$@" -d ;;
        --paired)     set -- "$@" -p ;;
        --help)       set -- "$@" -h ;;
        --version)    set -- "$@" -v ;;
        *)            set -- "$@" "$arg" ;;
    esac
done

# Wrapper function for help text
show_help() {
    echo "RNA-Seq pipeline helper, v${VERSION}"
    echo "Usage: $0 [options]"
    echo "  -f, --file         FILE  Specify input file with one SRA identifier per line"
    echo "  -d, --directory    DIRECTORY  Specify directory used for input/output. Must contain sub-directories 'Fastq' (with "
    echo "                     '.fastq.gz' files matching those given in <FILE>) and 'Genome' (with a '.fna' and '.gff3' file)"
    echo "  -p, --paired       Enable processing of paired-end data (with suffixes '_1' and '_2')"
    echo "  -h, --help         Show this help message and exit"
}

# If no arguments were provided, show help and exit
if [ $# -eq 0 ]; then
    show_help
    exit 0
fi

while getopts f:d:phv flag; do
    case "${flag}" in
        f) file=${OPTARG} ;;
        d) directory=${OPTARG} ;;
        p) p_flag=1 ;;
        h)
            show_help
            exit 0
            ;;
        v)
            echo ${VERSION}
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))


## Make sure the 'rename' command is installed
if ! command -v rename >/dev/null 2>&1; then
    echo "Command 'rename' not found. Check that 'rename' is installed"
    exit 1
fi


## Conda/Micromamba
source ${HOME}/micromamba/etc/profile.d/mamba.sh && micromamba activate rnaseq2025


gtf_file=$( ls ${directory}/Genome/*.gff3 )
fasta_file=$( ls ${directory}/Genome/*.fna )
dir_fastq=${directory}/Fastq
gtf_clean=${gtf_file%.gff3}_clean.gff3


## Clean the GFF3 file
sed "/##FASTA/q" ${gtf_file} | head -n -1 | grep -v "?" > ${gtf_clean}


## FastQC
## This will run quietly in the background by default; remove the trailing "&" to run it in the foreground, sequentially
if [ ! -d ${directory}/FastQC ]; then mkdir ${directory}/FastQC; fi
fastqc -t 4 -q -o ${directory}/FastQC/ ${directory}/Fastq/*.fastq.gz &
p1=$!


## Generate genome index
STAR \
    --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir ${directory}/Genome/ \
    --sjdbGTFfeatureExon CDS \
    --genomeFastaFiles ${fasta_file} \
    --sjdbGTFfile ${gtf_clean} \
    --genomeSAindexNbases 10


## STAR
if [ ! -d ${directory}/STAR ]; then mkdir ${directory}/STAR; fi

while read sample; do
  if [ "${p_flag:-0}" -eq 1 ]; then
    STAR \
      --runMode alignReads \
      --runThreadN 4 \
      --genomeDir ${directory}/Genome \
      --genomeLoad LoadAndKeep \
      --limitBAMsortRAM 48000000000 \
      --outSAMtype BAM SortedByCoordinate \
      --readFilesIn ${dir_fastq}/${sample}_1.fastq.gz ${dir_fastq}/${sample}_2.fastq.gz \
      --readFilesCommand zcat \
      --outFileNamePrefix ${directory}/STAR/${sample}_
  else
    STAR \
      --runMode alignReads \
      --runThreadN 4 \
      --genomeDir ${directory}/Genome \
      --genomeLoad LoadAndKeep \
      --limitBAMsortRAM 48000000000 \
      --outSAMtype BAM SortedByCoordinate \
      --readFilesIn ${dir_fastq}/${sample}.fastq.gz \
      --readFilesCommand zcat \
      --outFileNamePrefix ${directory}/STAR/${sample}_
  fi
done < ${file}

rename 's/_Aligned.sortedByCoord.out//' ${directory}/STAR/*.bam
STAR --genomeDir ${directory}/Genome --genomeLoad Remove
rm Aligned.out.sam Log.out Log.final.out Log.progress.out SJ.out.tab


## HTSeq
echo -en "\nStarting HTSeq..."
if [ ! -d ${directory}/HTSeq ]; then mkdir ${directory}/HTSeq; fi
find ${directory}/STAR/ -name "*.bam" | parallel --jobs 4 \
  "htseq-count -s reverse -a 10 -r pos -t CDS -i locus_tag {} "${gtf_clean}" > "${directory}"/HTSeq/{/.}.count"


## MultiQC
## Wait for FastQC to finish before running MultiQC
wait $p1
multiqc \
  -f \
  -o ${directory}/MultiQC \
  ${directory}/FastQC/ \
  ${directory}/STAR/ \
  ${directory}/HTSeq/


micromamba list > ${directory}/micromamba_env_rnaseq2025.txt
micromamba deactivate
