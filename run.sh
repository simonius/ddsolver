h5fc -march=native  -fopenmp -O3 consts.f90 miner.f90

ulimit -s unlimited
export OMP_STACKSIZE=4G
./a.out

