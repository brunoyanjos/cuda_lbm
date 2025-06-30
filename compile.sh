#!/bin/bash

# Prevent running with source
if [[ "$0" != "$BASH_SOURCE" ]]; then
    echo "Error: This script should not be run with 'source'. Please use:"
    echo "       bash compile.sh <simulation_id>"
    return 1 2>/dev/null || exit 1
fi

# Check if simulation ID is provided
if [ $# -eq 0 ]; then
    echo "Error: Simulation ID is required!" >&2
    echo "Usage: bash compile.sh <simulation_id>" >&2
    exit 1
fi

ID_SIM="$1"  # Required simulation ID
LT=D2Q9

# Set CompCap if not defined
if [ -z "$CompCap" ]; then
    CompCap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '.')
    if [ -z "$CompCap" ]; then
        echo "Error: Unable to determine compute capability." >&2
        exit 1
    fi
fi

# Clean specific executable
EXEC_NAME="sim_${LT}_sm${CompCap}"
rm -f "$EXEC_NAME"

# Compile from source directory
cd MR_LBM/src || {
    echo "Error: Failed to enter MR_LBM/src directory" >&2
    exit 1
}

nvcc -gencode arch=compute_${CompCap},code=sm_${CompCap} -rdc=true -O3 --restrict \
    -D ID_SIM=\"$ID_SIM\" \
     *.cu \
     -lcudadevrt -lcurand -o ../../"$EXEC_NAME" || {
    echo "Error: Compilation failed" >&2
    cd ../..
    exit 1
}

# Return and run with simulation ID
cd ../..
./"$EXEC_NAME" || {
    echo "Error: Simulation failed" >&2
    exit 1
}