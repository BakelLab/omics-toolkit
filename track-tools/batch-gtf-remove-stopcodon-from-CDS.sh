#!/bin/bash

# 14.06.2020 11:55:14 EDT
# Harm van Bakel <hvbakel@gmail.com>

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail

# --- Check usage ---
if (( $# != 2 )); then
    echo -e "\n  Usage: $(basename "$0") <isoseq.gtf> <output.gtf>\n" >&2
    exit 1
fi

# --- Create a secure temporary directory ---
# The trap command will ensure this directory is cleaned up on EXIT.
# The 'mktemp -d' command creates a unique directory, avoiding name collisions.
TEMP_DIR=$(mktemp -d)

# --- Setup a trap to clean up the temporary directory on exit ---
# This function will run when the script exits, for any reason (success, error, or interrupt).
cleanup() {
    # The '|| true' prevents the script from erroring out if the directory is already gone.
    rm -rf "${TEMP_DIR}" || true
}
trap cleanup EXIT

# --- Split input file by chromosome ---
# Using 'mapfile' (or 'readarray') is a safe way to read lines into an array.
mapfile -t chr_list < <(cut -f 1 "$1" | grep -vP '^\s*#' | sort -u)

echo "Processing ${#chr_list[@]} chromosomes..."

for i in "${chr_list[@]}"; do
    awk -F '\t' -v OFS='\t' -v chr="$i" '$1==chr' "$1" > "${TEMP_DIR}/input_${i}.gtf"
done

# --- Run processing in parallel ---
# All temporary files are now safely contained within TEMP_DIR.
# We pass the array of chromosomes to parallel using "${chr_list[@]}".
parallel --jobs 24 "gtf-remove-stopcodon-from-CDS.R -i ${TEMP_DIR}/input_{}.gtf -o ${TEMP_DIR}/output_{}.gtf 2> ${TEMP_DIR}/output_{}.log" ::: "${chr_list[@]}"

# --- Combine final outputs ---
# The glob will expand to find all output files inside the temporary directory.
cat "${TEMP_DIR}"/output_*.gtf > "$2"

echo "Processing complete. Output written to $2"

# The 'trap cleanup EXIT' command will automatically handle cleanup from here.
