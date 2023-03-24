#!/bin/bash
#PJM -L "node=1"
#PJM -L "rscgrp=small"
#PJM -L "freq=2000"
#PJM -L "elapse=3:30:00"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -S

. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load /zc5pwgc

mkdir result

algo=("seq" "gnu" "psrs" "pses")

p=7
N=$((10**p))
count=20

# "sequential algo name", "parallel algo number", "num threads"
test(){
    i=$2
    export OMP_NUM_THREADS=$3
    fipp -C -d profile -Icall,cpupa -m 10000 ./${1}.out ${N} ${count} ${i}>> "result/${algo[i]}_${1}${p}.time" 2>> "result/${algo[i]}_${1}${p}.err"
    fipppx -A -pall -Ibalance,call -d profile -o "result/${algo[i]}_${1}${p}_${3}.prof"
    rm -rf profile
}

test 'blockq' 0 1
test 'std' 0 1
test 'pdq' 0 1

for th in 1 2 4 12 24 48
do
    export XOS_MMM_L_PAGING_POLICY=demand:demand:demand
    numactl -C 12-$((11+th)) -m 4-$((4+(th-1)/12))

    test 'blockq' 1 ${th}

    for i in 2 3
    do
        test 'blockq' ${i} ${th}
        test 'std' ${i} ${th}
        test 'pdq' ${i} ${th}
        test 'std_heap' ${i} ${th}
        test 'std_sort' ${i} ${th}
    done
done
