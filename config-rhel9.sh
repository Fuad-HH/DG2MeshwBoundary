#SIM_VER=simmetrix-simmodsuite-2024.1-240620dev-g4avy72vpgsyvhc4wkwbg4zh2kn4v74z
#SIM_VER=simmetrix-simmodsuite-2025.0-241016dev-vafjs2qz5n2rokmrhs5ijvsua64ty4kz
SIM_VER=simmetrix-simmodsuite-2025.1-250602dev-yv5oiomk2sdcru5cn5vqjxprycytc4ow
SIM_ARCHOS=x64_rhel8_gcc83

cmake -S . -B build \
	-DSIM_MPI=mpich4.1.1 \
	-DSIMMETRIX_INCLUDE_DIR=/opt/scorec/spack/rhel9/v0201_4/install/linux-rhel9-x86_64/gcc-12.3.0/$SIM_VER/include \
	-DSIMMETRIX_LIB_DIR=/opt/scorec/spack/rhel9/v0201_4/install/linux-rhel9-x86_64/gcc-12.3.0/$SIM_VER/lib/$SIM_ARCHOS \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DSIM_PARASOLID=ON

cmake --build build -j2
