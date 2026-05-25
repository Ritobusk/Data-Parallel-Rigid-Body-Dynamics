#!/bin/bash
#SBATCH --job-name=rnea
# normal cpu stuff: allocate cpus, memory
#SBATCH --output=dynamics_bench_a100_cpu_and_blocked.out
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=80000M
#SBATCH -p gpu --gres=gpu:a100:1
#SBATCH --time=2:00:00

#your script, in this case: write the hostname and the ids of the chosen gpus.
lscpu
nvidia-smi
nvcc --version
hostname
echo $CUDA_VISIBLE_DEVICES
echo "\nRNEA vtree GPU - Optimal Data Structure\n"
futhark bench --backend=cuda bench_rnea_optimal_ds.fut
echo "\nRNEA vtree GPU\n"
futhark bench --backend=cuda bench_rnea2.fut
echo "\nRNEA CPU\n"
futhark bench --backend=c bench_rnea_no_vtree.fut

echo "\nRNEA vtree GPU - Optimal Data Structure\n"
futhark bench --backend=cuda bench_crba_optimal_ds.fut
echo "\nCRBA vtree GPU\n"
futhark bench --backend=cuda bench_crba2.fut
echo "\nCRBA CPU\n"
futhark bench --backend=c bench_crba2_seq.fut
