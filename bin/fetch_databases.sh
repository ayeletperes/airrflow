#!/usr/bin/env bash
# Written by Gisela Gabernet and released under the MIT license (2020).
# Reworked to fetch and build with sourcerer instead of the fetch/ref2igblast
# shell scripts. sourcerer's source names match the database types (imgt,
# airrc-imgt), so the type is passed through directly.
set -euo pipefail

DATABASE_TYPE=imgt

while getopts "d:" OPT; do
    case "$OPT" in
    d)  DATABASE_TYPE=$OPTARG
        ;;
    esac
done

echo "Fetching databases with sourcerer (${DATABASE_TYPE})..."

# Download each species into one reference tree. sourcerer writes the germline
# FASTAs under <outdir>/reference_base/<species>/, so both species accumulate in
# the same reference_base.
for SPECIES in human mouse; do
    sourcerer "${DATABASE_TYPE}" download "${SPECIES}" --outdir sourcerer_out
done
mv sourcerer_out/reference_base reference_base

# Build the IgBLAST databases from the reference tree. This cleans the FASTAs,
# runs makeblastdb, and mirrors the NCBI internal_data / optional_file trees into
# igblast_base -- the work fetch_igblastdb.sh and ref2igblast.sh used to do.
sourcerer reference build reference_base --out igblast_base

echo "FetchDBs process finished."
