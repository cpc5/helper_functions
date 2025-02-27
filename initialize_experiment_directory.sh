#!/bin/bash

# Parse command line arguments
while getopts "n:" opt; do
    case $opt in
        n) dirname="$OPTARG";;
        *) echo "Usage: $0 -n <path/experimentname>" >&2
           exit 1;;
    esac
done

# Check if dirname was provided
if [ -z "$dirname" ]; then
    echo "Error: Directory name must be provided with -n option"
    echo "Usage: $0 -n <path/experimentname>"
    exit 1
fi

# Create main directory
mkdir -p "$dirname"

# Create subdirectories
mkdir -p "$dirname/code"
mkdir -p "$dirname/data"
mkdir -p "$dirname/pdfs"

chmod 750 "$dirname" -R
chgrp pi-imoskowitz "$dirname" -R
echo "Created directory structure in $dirname"
