# DG2MeshwBoundary

Creates unstructured Degas2 mesh with the given boundary nodes and edges to conform XGC mesh boundary.
---
## Build Instructions
In SCOREC RHEL9 machines, load the following modules:
```bash
module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
module load gcc/12.3.0-iil3lno mpich/4.1.1-xpoyz4t
module load simmetrix-simmodsuite/2025.1-250602dev-yv5oiom
```
Then configure and build using CMake:
```bash
cmake -S . -B build \
        -DSIM_MPI=mpich4.1.1 \
        -DSIMMETRIX_INCLUDE_DIR=$SIMMETRIX_SIMMODSUITE_ROOT/include/ \
        -DCMAKE_CXX_COMPILER=mpicxx \
        -DCMAKE_C_COMPILER=mpicc \

cmake --build build -j2
```

or just use the provided build script:
```bash
. config-rhel9.sh
```

## Usage
See examples in [EXAMPLES.md](EXAMPLES.md).