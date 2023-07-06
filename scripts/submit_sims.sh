#!/bin/bash
#Submit assembly simulations

#SBATCH --constrain="lustre"
#SBATCH --account=TG-MCB090163
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH -N 1 
#SBATCH -n 1

in_file=$1
seed=$2
seed_file=$3

echo "Running assembly simulation with input file: $in_file"

run_dir="/home/lfrechette/lj_assembly/bin/"

$run_dir/lj_assembly $in_file $seed $seed_file

