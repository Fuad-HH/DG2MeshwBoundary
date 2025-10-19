cmake -S . -B build \
	-DSIM_MPI=mpich4.1.1 \
	-DSIMMETRIX_INCLUDE_DIR=$SIMMETRIX_SIMMODSUITE_ROOT/include/ \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \

cmake --build build -j2
