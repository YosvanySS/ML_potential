# xyz2qe
The main goal of this Python script is to split the database for any Machine Learning Potential (MLP) in the format `*.xyz` into their respective input files of [Quantum Espresso (QE)](https://www.quantum-espresso.org/). The script executes QE locally for 2 seconds to determine the number of _k-points_ and Fast Fourier Transform coefficients in the Z direction (FFTz), which are used for parallelizing the code in the Cluster.

> [!IMPORTANT]  
> <ins>**Requirements:**</ins>\
> `QE` installed. See [home page](https://www.quantum-espresso.org/) for downloading it.\
> `timeout` installed. It kills QE's execution after a 2s in our local machine.\
> `parallel` installed. It runs QE simultaneously on different CPUs in our local machine.

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

