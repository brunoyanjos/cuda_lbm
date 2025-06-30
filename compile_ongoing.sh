#!/bin/bash

# Prevent running with source
if [[ "$0" != "$BASH_SOURCE" ]]; then
    echo "Error: This script should not be run with 'source'. Please use:"
    echo "    bash compile.sh"
    echo "    bash compile.sh <simulation_id> <case_type>"
    return 1 2>/dev/null || exit 1
fi

# Check arguments count
if [ $# -eq 2 ]; then
    # Non-interactive mode
    SIM_ID="$1"
    CASE_TYPE="$2"
elif [ $# -eq 0 ]; then
    # Interactive mode
    echo "┌──────────────────────────────────────────────┐"
    echo "│          Simulation Setup Wizard             │"
    echo "└──────────────────────────────────────────────┘"
    
    # Get simulation ID
    while true; do
        echo -n "Enter simulation ID (e.g., 001): "
        read SIM_ID
        if [ -n "$SIM_ID" ]; then
            break
        else
            echo "Error: ID cannot be empty!"
        fi
    done
    
    # Get case type
    while true; do
        echo "Select case type:"
        echo "1) CYLINDER"
        echo "2) LDC (Lid-Driven Cavity)"
        echo "3) JET_FLOW"
        echo -n "Enter choice (1-3): "
        read CHOICE
        
        case $CHOICE in
            1) CASE_TYPE="CYLINDER"; break;;
            2) CASE_TYPE="LDC"; break;;
            3) CASE_TYPE="JET_FLOW"; break;;
            *) echo "Invalid choice! Please enter 1-3.";;
        esac
    done
else
    echo "Error: Invalid number of arguments!" >&2
    echo "Usage:" >&2
    echo "  Interactive mode: bash compile.sh" >&2
    echo "  Non-interactive:  bash compile.sh <ID> <CASE_TYPE>" >&2
    echo "Valid CASE_TYPE options: CYLINDER, LDC, JET_FLOW, BACKSTEP, POROUS_MEDIUM" >&2
    exit 1
fi

# Validate case type
VALID_CASES=("CYLINDER" "LDC" "JET_FLOW")
if [[ ! " ${VALID_CASES[@]} " =~ " ${CASE_TYPE} " ]]; then
    echo "Error: Invalid case type '$CASE_TYPE'!" >&2
    echo "Valid options: ${VALID_CASES[*]}" >&2
    exit 1
fi

echo "──────────────────────────────────────────────"
echo " Starting $CASE_TYPE simulation with ID: $SIM_ID"
echo "──────────────────────────────────────────────"

# Set Lattice Type
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

# Pass parameters to compiler
nvcc -gencode arch=compute_${CompCap},code=sm_${CompCap} -rdc=true -O3 --restrict \
     -D CASE_TYPE=\"$CASE_TYPE\" \
     *.cu \
     -lcudadevrt -lcurand -o ../../"$EXEC_NAME" || {
    echo "Error: Compilation failed" >&2
    cd ../..
    exit 1
}

# Return and run with simulation parameters
cd ../..
./"$EXEC_NAME" "$SIM_ID" "$CASE_TYPE" || {
    echo "Error: Simulation failed" >&2
    exit 1
}