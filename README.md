# Julia framework for automated simulations in COM3

This framework was developed for the purpose of easily simulating compression tests of PVA-ECC-based 2D metamaterials using COM3.\
It contains all the necessary tools to generate truss-like geometries, mesh them, write COM3 files, run instances of COM3 and read from COM3 output files,
as well as post-processing of the data such as exporting results to `.csv` and producing plots. It also contains a Reinforcement Learning environment built upon `Reinforce.jl`
which allows for easy implementation of reinforcement learning algorithms using COM3 simulations to calculate rewards.\
The whole framework is made thread-safe so that many simulations can be run in parallel on CPUs with many cores, improving performance on big simulation batches.

## Starting up with the framework

This framework is a module coded using the *Julia programming language*, a high performance language for scientific purposes.
You will need to first install *Julia* (https://julialang.org/) and *Git for Windows* (https://gitforwindows.org/). Once that is done, use the *git GUI* to clone this github repository.\
Go to the windows command line, ( by typing `cmd` in the search bar for instance) and type in :
```
setx PATH "%PATH%;C:\Users\[Your-User-Name]\AppData\Local\Programs\Julia-1.6.1\bin"
```
after replacing with your correct username.\
Then launch *git bash* and run the following in the Metamaterials folder (use `cd` to change directories inside the bash, `pwd` to display the current directory) :
```
julia setup.jl
``` 
This script will install relevant packages for this framework.\
Once this is done, you can run a simple batch of simulations by writing in the *git bash* terminal, making sure you are in the `Metamaterials` directory :
```
julia batchSimulator.jl -batchsize N
```
where `N` is the number of simulations you want to run.\
For an overview of the supported arguments by the framework, type :
```
julia batchSimulator.jl -help
```
If you want to run simulations simultameously across several threads, type :
```
julia --threads (Your desired number of threads) batchSimulator.jl (Further options...)
```
Note that this is generally faster as the simulations are quite small scale, thus not taking advantage of many processing threads. Running instances in parallel, however, allows for faster computation.

## Structure of the Julia module

The module is comprised of several source files, each defining functions for a specific part of the framework :
* `structures.jl` defines the data structures used in the module (*Point,Element,Skeleton,Model,Material,Simulation*) along with various constructors an methods for them.
* `geometryTools.jl` contains the geometry generator as well as plotting functions to this end.
* `simulationTools.jl` holds the tools which build and run a simulation.
* `COM3RWTools` gives the different methods for writing to and reading from COM3 files.
