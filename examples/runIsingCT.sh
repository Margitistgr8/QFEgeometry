#!/bin/bash

EXEC=/Users/jinyunlin/QFE-Research/QFEgeometry/bin/IsingflatCT

Nx=40
Ny=40
invT=5.0
k1=1.0
k2=1.0
n_iter=50000
n_therm=5000
n_skip=10
baseseed=12342
i=0
#couplings=(2.95 2.975 3.025)
couplings=(1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5)

# Parallelism control
export OMP_NUM_THREADS=2   # each simulation uses 2 threads
MAX_JOBS=3                 # max number of concurrent simulations

# Internal counters
i=0
job_count=0

OUTPUT_DIR="/Users/jinyunlin/QFE-Research/QFEgeometry/data/IsingCT"
# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

for tc in "${couplings[@]}"; do
    seed=$((baseseed + i))
    echo "Launching simulation for tc=$tc with seed=$seed..."

    # Run simulation with timeout and per-job log file
    timeout 3h $EXEC -x $Nx -y $Ny \
          -T $invT -a $k1 -b $k2 \
          -c $tc -h $n_iter -t $n_therm -s $n_skip \
          -S $seed -d $OUTPUT_DIR > "${OUTPUT_DIR}/log_tc_${tc}_${Nx}.out" 2>&1 &

    ((job_count++))
    ((i++))

    # If max parallel jobs reached, wait for all to finish before continuing
    if (( job_count >= MAX_JOBS )); then
        wait
        job_count=0
    fi
done

# Wait for any remaining background jobs to finish
wait

echo "✅ All simulations completed."

# for tc in 0.05 0.1 0.2 0.4 0.8
# for tc in 2.9
# do
#     seed=$((baseseed+i))
#     $EXEC -x $Nx -y $Ny \
#           -T $invT -a $k1 -b $k2 \
#           -c $tc -h $n_iter -t $n_therm -s $n_skip \
#           -S $seed -d $OUTPUT_DIR
# done