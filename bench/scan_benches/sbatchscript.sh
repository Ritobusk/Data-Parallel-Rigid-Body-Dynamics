#!/bin/bash
#SBATCH --job-name=sc_variations
# normal cpu stuff: allocate cpus, memory
#SBATCH --output=vtee_scan_benches_a100.out
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=40000M
#SBATCH -p gpu --gres=gpu:a100:1
#SBATCH --time=1:00:00

#your script, in this case: write the hostname and the ids of the chosen gpus.
nvidia-smi
#nvcc --versin
hostname
echo $CUDA_VISIBLE_DEVICES
futhark bench --backend=cuda scan_variations.fut
futhark bench --backend=cuda we_scan_bench.fut
futhark bench --backend=cuda blocked_scan_bench.fut
