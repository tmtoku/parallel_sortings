#!/bin/bash

# boost library
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load /zc5pwgc

cur_dir=$(pwd)

# std.out
cd ../src &&
make clean && seq=std make &&
cp a.out ${cur_dir}/std.out &&
cd ${cur_dir}

# pdq.out
cd ../src &&
make clean && seq=pdq make &&
cp a.out ${cur_dir}/pdq.out &&
cd ${cur_dir}

# blockq.out
cd ../src &&
make clean && make &&
cp a.out ${cur_dir}/blockq.out &&
cd ${cur_dir}

# std_heap.out
cd ../src &&
make clean && seq=std merge=heap make &&
cp a.out ${cur_dir}/std_heap.out &&
cd ${cur_dir}

# std_sort.out
cd ../src &&
make clean && seq=std merge=sort make &&
cp a.out ${cur_dir}/std_sort.out &&
cd ${cur_dir}
