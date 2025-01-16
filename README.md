# xyz2qe
The main goal of this Python script is to split the database for any Machine Learning Potential (MLP) in the format `*.xyz` into their respective input files of [Quantum Espresso](https://www.quantum-espresso.org/) (QE).\
The script executes QE locally for 2s to determine the number of _k-points_ and other considered parameters for the paralellization of the code in the Cluster.

> [!IMPORTANT]  
> **Requirements for using this script:**\
> `timeout` installed. It is used to killed Quantum Expresso (QE) after a while running in our local machine.\
> `parallel` installed. It is used to run QE at the same time in different CPUs in our local machine.

> [!TIP]
> 
> **Installing `timeout` and `parallel` in MacOS via MacPorts**
> ```
> sudo port install timeout
> sudo port install parallel
> ```
>
> **Installing `timeout` and `parallel` in Linux**
> ```
> sudo apt-get install timeout
> sudo apt-get install parallel
> ```

> [!NOTE]
> The main goal of this Python script is to split the database of any file in the format `*.xyz` into their respective input files of QE.\
> The script has to be placed in the same folder of the database file.
>
> 
> 
> 'input_qe' is reading the database and creating input files for Quantum Espresso (QE) per each structure in the database.
> It also creates a temporary bash repository to run QE in the local machine to determine the parameters of parallelization.
> It creates two folders: 'local' and 'cluster'. In 'local' is run QE while in 'cluster' are placed the run files and inputs for running in the cluster.

> [!WARNING]  
> Critical content demanding immediate user attention due to potential risks.

> [!CAUTION]
> Negative potential consequences of an action.

