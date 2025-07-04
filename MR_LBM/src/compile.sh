LT=D2Q9

# Set CompCap if it's not already defined
if [ -z "$CompCap" ]; then
    CompCap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '.')
    if [ -z "$CompCap" ]; then
        echo "Error: Unable to determine compute capability."
        exit 1
    fi
fi
rm -f ./../*sim_D2Q9_sm75
# rm -f ./2D_CHANNEL/$1/*.dat
# rm -r ./2D_CHANNEL/$1
    nvcc -gencode arch=compute_${CompCap},code=sm_${CompCap} -rdc=true --ptxas-options=-v -O3 --restrict \
        *.cu \
        -lcudadevrt -lcurand -o ./../$1sim_${LT}_sm${CompCap}
                
./../$1sim_${LT}_sm${CompCap}     



