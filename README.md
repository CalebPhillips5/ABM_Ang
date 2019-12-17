# ABM_Ang
Mathematical model of tumor angiogenesis
Quite a few libraries will be required to utilize the code "ABM_unbreakable"
Petsc, Libmesh (and other various valgrind gnu mpich gsl boost)

The smaller tools generally come with Xcode on Mac or are built in on other platforms.

Download Petsc source files
cd /source/petsc-3.8.2/

% Configure petsc
./configure --prefix=/PETSc382/ --download-cmake --with-mpi-dir=$MPI_DIR --with-debugging=false --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3' --download-fblaslapack=1 --with-mumps=1 --download-mumps=1 --with-metis=1 --download-metis=1 --with-parmetis=1 --download-parmetis=1 --with-scalapack=1 --download-scalapack=1 --with-superlu=1 --download-superlu=1 --with-valgrind-dir=/opt/ohpc/pub/valgrind/3.11.0/
make PETSC_DIR=/data/source/petsc-3.8.2 PETSC_ARCH=arch-linux2-c-opt all
make PETSC_DIR=/data/source/petsc-3.8.2 PETSC_ARCH=arch-linux2-c-opt install

Download Libmesh source files
cd /source/libmesh-1.2.1

% Configre Libmesh
export PETSC_DIR=/PETSc382/
./configure --prefix=/libMesh121/ --with-metis=PETSc
make
make install

%% Running the Code!
full_abm.c        - this allows you to change some of the parameters we studied in a sensitivity analysis
main_abm_oxy.c    - the main function that execudes the code
cell_abm.c        - houses both update_states and compute_forces, the crux of the simulations
  update_state    - where phenotypic transitions take place
  compute_forces  - where the forces are computed and the cells move accordingly
generate_figures  - generates the ABM figures
options.in        - simulation settings (how many timesteps, where to restart, print solution, etc.)

