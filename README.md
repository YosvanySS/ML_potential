# ML_potential

> [!IMPORTANT]  
> Requirements for using this script:
> `timeout` installed. It is used to killed Quantum Expresso (QE) after a while running in our local machine.
> `parallel` installed. It is used to run QE at the same time in different CPUs in our local machine.
>
> [!TIP]
> ** Installing `timeout` and `parallel` in MacOS via MacPorts **
> ```
> sudo port install timeout
> sudo port install parallel
> ```

> [!WARNING]  
> Critical content demanding immediate user attention due to potential risks.

> [!CAUTION]
> Negative potential consequences of an action.


> [!NOTE]
> 'input_qe' is reading the database and creating input files for Quantum Espresso (QE) per each structure in the database.
> It also creates a temporary bash repository to run QE in the local machine to determine the parameters of parallelization.
> It creates two folders: 'local' and 'cluster'. In 'local' is run QE while in 'cluster' are placed the run files and inputs for running in the cluster.
