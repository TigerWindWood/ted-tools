#!/bin/bash

# This file is a part of TED: The Encyclopedia of Domains. If you utilize or reference any content from this file, 
# please cite the following paper:

# Lau et al., 2024. Exploring structural diversity across the protein universe with The Encyclopedia of Domains.

# Function to display usage message
usage() {
    echo "Usage: $0 -i <input_directory_with_pdb_files> -o <output_directory>"
    exit 1
}

# Check that the environment exists and activate it 
BASE_DIR="/root/data/ted-tools-main/ted_consensus_1.0"
VENV_DIR="${BASE_DIR}/ted_consensus"
if [ -d "$VENV_DIR" ]; then
    source $VENV_DIR/bin/activate
else
    echo "Virtual environment 'ted_consensus' does not exist."
    echo "Please run 'bash setup.sh' to create and set up the virtual environment."
    exit 1
fi

# Parse command-line arguments
while getopts "i:o:" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if both input and output directories are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

# Check if the input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: $INPUT_DIR is not a directory"
    exit 1
fi

# Create the output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

SCRIPT_DIR="/root/data/ted-tools-main/ted_consensus_1.0"
PY=$(which python)

SEGMENT="${SCRIPT_DIR}/scripts/segment.sh"
CONSENSUS="${SCRIPT_DIR}/scripts/get_consensus.py"
FILTER_DOMAINS="${SCRIPT_DIR}/scripts/filter_domains_consensus.py"

# Calculate input data
input_pdb_count=$(find "${INPUT_DIR}" -maxdepth 1 -name '*.pdb' | wc -l)

# Run Merizo on the input directory
out_merizo="${OUTPUT_DIR}/chopping_merizo.txt"
log_merizo="${OUTPUT_DIR}/chopping_merizo.log"
if test -f "${out_merizo}"; then
    # Output file exists, check the line count
    merizo_count=$(wc -l < "${out_merizo}")
    # Result count is equal, skip execution
    if [ "${merizo_count}" -eq "${input_pdb_count}" ]; then
        echo "${out_merizo} already exists"
    # Result count is not equal, delete the output file and log file, and execute again
    else
        rm -f "${out_merizo}" "${log_merizo}"
        bash "${SEGMENT}" -i "${INPUT_DIR}" -m merizo -o "${OUTPUT_DIR}" > "${log_merizo}" 2>&1
    fi
# Output file does not exist, execute directly
else
    rm -f "${out_merizo}" "${log_merizo}"
    bash "${SEGMENT}" -i "${INPUT_DIR}" -m merizo -o "${OUTPUT_DIR}" > "${log_merizo}" 2>&1
fi

if test ! -f "${out_merizo}" || test ! -s "${out_merizo}"; then
    echo "Expected to find chopping file for Merizo at ${out_merizo}!"
    exit 1
fi

# Run UniDoc on the Merizo output
out_unidoc="${OUTPUT_DIR}/chopping_unidoc.txt"
log_unidoc="${OUTPUT_DIR}/chopping_unidoc.log"
if test -f "${out_unidoc}"; then
    # Output file exists, check the line count
    unidoc_count=$(wc -l < "${out_unidoc}")
    if [ "${unidoc_count}" -eq "${input_pdb_count}" ]; then
        echo "${out_unidoc} already exists"
    else
        rm -f "${out_unidoc}" "${log_unidoc}"
	bash "${SEGMENT}" -i "${INPUT_DIR}" -m unidoc -o "${OUTPUT_DIR}" -c "${out_merizo}" > "${log_unidoc}" 2>&1
    fi
else
    rm -f "${out_unidoc}" "${log_unidoc}"
    bash "${SEGMENT}" -i "${INPUT_DIR}" -m unidoc -o "${OUTPUT_DIR}" -c "${out_merizo}" > "${log_unidoc}" 2>&1
fi

if test ! -f "${out_unidoc}" || test ! -s "${out_unidoc}"; then
    echo "Expected to find chopping file for UniDoc at ${out_unidoc}!"
    exit 1
fi

# Run Chainsaw on the input directory
out_chainsaw="${OUTPUT_DIR}/chopping_chainsaw.txt"
log_chainsaw="${OUTPUT_DIR}/chopping_chainsaw.log"
if test -f "${out_chainsaw}"; then
    # Output file exists, check the line count
    chainsaw_count=$(wc -l < "${out_chainsaw}")
    if [ "${chainsaw_count}" -eq "${input_pdb_count}" ]; then
        echo "${out_chainsaw} already exists"
    else
        rm -f "${out_chainsaw}" "${log_chainsaw}"
        bash "${SEGMENT}" -i "${INPUT_DIR}" -m chainsaw -o "${OUTPUT_DIR}" > "${log_chainsaw}" 2>&1
	fi
else
    rm -f "${out_chainsaw}" "${log_chainsaw}"
    bash "${SEGMENT}" -i "${INPUT_DIR}" -m chainsaw -o "${OUTPUT_DIR}" > "${log_chainsaw}" 2>&1
fi

if test ! -f "${out_chainsaw}" || test ! -s "${out_chainsaw}"; then
    echo "Expected to find chopping file for Chainsaw at ${out_chainsaw}!"
    exit 1
fi


echo "Calculating consensus domains from Merizo, UniDoc and Chainsaw outputs.. "
# Calculate consensus from each of the outputs
out_consensus="${OUTPUT_DIR}/consensus.tsv"
log_consensus="${OUTPUT_DIR}/consensus.log"
"${PY}" "${CONSENSUS}" -c "${out_merizo}" "${out_chainsaw}" "${out_unidoc}" -o "${out_consensus}" > "${log_consensus}" 2>&1

if test -f "${out_consensus}"; then
    "${PY}" "${FILTER_DOMAINS}" "${out_consensus}" -o "${out_consensus}.tmp"

    if [ $? == 0 ]; then
        mv "${out_consensus}.tmp" "${out_consensus}"
    fi
else
    echo "Expected to find consensus domain file at ${out_consensus}"
    exit 1
fi

echo "Consensus domain file saved at ${out_consensus}"
