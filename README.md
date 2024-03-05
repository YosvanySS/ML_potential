# ML_potential

> [!NOTE]
> 'input_qe' is reading the database and creating input files for Quantum Espresso (QE) per each structure in the database.
> It also creates a temporary bash repository to run QE in the local machine to determine the parameters of parallelization.
> It creates two folders: 'local' and 'cluster'. In 'local' is run QE while in 'cluster' are placed the run files and inputs for running in the cluster.
