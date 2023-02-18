#!/bin/bash

METRICA="L3 L2CACHE FLOPS_DP FLOPS_AVX"
TAMANHOS="34 32 66"

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

for n in $TAMANHOS
do
    for k in $METRICA
    do
        likwid-perfctr -C 3 -g ${k} -m -O ./cgSolver -n ${n} -k 7 -p 1 -i 150 -o outputs/output_${n}_${k}.txt > logs/${k}_${n}.log
    done
done

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor 
