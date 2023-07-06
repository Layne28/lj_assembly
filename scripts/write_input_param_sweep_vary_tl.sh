#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
va=$5 
kT=$6 

dt=0.0002
freq=1000
equil_steps=50000
production_steps=5000000

base_out_dir="/expanse/lustre/projects/brs104/lfrechette/"

taus=(0.01 0.10 0.30 1.00 10.00)
lambdas=(0.50 1.00 2.00 3.00 5.00 10.00)
phis=(0.20 0.40 0.60)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for phi in "${phis[@]}"
        do
            python3 scripts/write_input_file.py input_files --va $va --tau $tau --Lambda $lambda --kT $kT --phi $phi --Lx $Lx --Ly $Ly --nx $nx --ny $ny --particles_freq $freq --thermo_freq $freq --dt $dt --equil_steps=$equil_steps --production_steps=$production_steps --particle_protocol "lattice" --output_dir "${base_out_dir}/lj-assembly-results/kT=${kT}/phi=${phi}/va=${va}/tau=${tau}/lambda=${lambda}/Lx=${Lx}_Ly=${Ly}/nx=${nx}_ny=${ny}"
        done
    done
done
