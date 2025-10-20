#!/bin/bash
set -e
shopt -s nullglob # To prevent errors in the for loop when there are no files left
VERSION="0.1.0"

## You MUST pass the flags for these arguments, e.g. "./helper_download_sra.sh -f samples.txt -d GSE151953"
for arg in "$@"; do
    shift
    case "$arg" in
        --file)      set -- "$@" -f ;;
        --directory) set -- "$@" -d ;;
        --help)      set -- "$@" -h ;;
        --version)   set -- "$@" -v ;;
        *)           set -- "$@" "$arg" ;;
    esac
done

# Help information wrapped into a function
show_help() {
    echo "SRA download helper, v${VERSION}"
    echo "Usage: $0 [options]"
    echo "  -f, --file FILE            Input file with one SRA identifier per line"
    echo "  -d, --directory DIRECTORY  Directory used for input/output"
    echo "  -h, --help                 Show this help message and exit"
    echo "  -v, --version              Show version information and exit"
}

# If no arguments were provided, show help and exit
if [ $# -eq 0 ]; then
    show_help
    exit 0
fi

while getopts f:d:hv flag; do
    case "${flag}" in
        f) file=${OPTARG};;
        d) directory=${OPTARG};;
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


## Veryify needed programs are installed
if ! command -v prefetch > /dev/null 2>&1; then
    echo "Command 'prefetch' not found. Check that 'sra-toolkit' is installed"
    exit 1
fi

if ! command -v fasterq-dump > /dev/null 2>&1; then
    echo "Command 'fasterq-dump' not found. Check that 'sra-toolkit' is installed"
    exit 1
fi

if ! command -v pigz > /dev/null 2>&1; then
    echo "Command 'pigz' not found. Check that 'pigz' is installed"
    exit 1
fi


echo "Downloading and compressing SRA data..."
process_acc() {
    acc="$1"
    prefetch ${acc} && fasterq-dump -m 4G -p ${acc}
    pigz -p 4 ${acc}*.fastq
    rm -r ${acc}
}
export -f process_acc
cat ${file} | parallel --verbose -j 4 process_acc

echo "Cleaning up..."
mkdir -p ${directory}/Fastq
mv *.fastq.gz ${directory}/Fastq/

echo "Done."
