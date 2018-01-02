# Energy Infrastructure Adaptive Design for Climate Change
[![Build Status](https://ci.lanlytics.com/sitew/ADCC.svg?token=TPN7vp19bq64KY6rQeyS&branch=master)](https://ci.lanlytics.com/sitew/ADCC)
[![codecov](https://cov.lanlytics.com/ghe/ansi/PowerModelsLANL.jl/branch/master/graph/badge.svg)](https://cov.lanlytics.com/ghe/sitew/ADCC)

This project is designed for exploring power systems incremental designs in response to future climate changes. Potential global climate change is happening and could potentially bring gradual/drastic climate changes, such as sea level raises, frequent storm surges, etc. These changes brings challenges upon current power systems designs, especially in coastal regions. Besides, such climate changes can also re-shape the current established power profiles? which indirectly affects the future design of the power systems. The sophistication of designing future grids incorporating climate change considerations lays in the uncertainty of these changes as well as the climate consequences they brought.

## Installation
Clone this repo to your local environment.
```git clone https://github.lanlytics.com/sitew/ADCC```

It is important to deploy an configuration file (.json) to control the algorith without much overhead on command line control.
You will find an example `config.json`.
The algorithm search both home dir
```>$HOME/Debug/Climate/config.json```and repo dir```>./config.json```
for this config.json

## Dependencies
This small tool requires the following packages in julia

```julia
> Pkg.add("JuMP")
> Pkg.add("JSON")
> Pkg.add("PowerModels")
> Pkg.add("ProgressMeter")
> Pkg.add("Distributions")
> Pkg.add("ArgParse")
> Pkg.add("Compat")
> Pkg.add("Glob")
```

Additional packages might be required in case what MIP solver is used
```julia
> Pkg.add("Gurobi") # If using it
> Pkg.add("CPLEX") # If using it
```

## To do a adaptive design optimization
Let's try this simple command with example scenario inputs. When you are in the repo,
```julia
> julia adcc.jl ieee118 dc regular --sm=file --sf="./format/scenario_example.json" --T=5 --EPS=0.0 --pf=paper_param --SOLVER=Cplex --TIMELIMIT=300 --NAME=yourscenario
```
The above command will conduct an adaptive design optimization given IEEE-118 (`ieee118`) test network
with DC dispatch mode (`dc`) using directly CPLEX (`regular` & `--SOLVER=Cplex`) with 300s time limits (`--TIMELIMIT=300`).
Note that `ieee118`, `dc`, and `regular` are positional inputs to specify a network, a dispatch mode, and a algorithm for optimization.
The optimization is conducted with the scenario input file (`--sm=file`) stored within the repo (`--sf=./format/scenario_example.json`).
Argument `--EPS=0.0` ensures that no risk is allowed, which means all scenarios in `scenario_example.json` should be satisfied.
The output of this run will be named `eval_yourscenario.json` (`--NAME=yourscenario`), where prefix and postfix are all added.
This output file will be find through `OUTPUTPATH` in your configuration file. Similarly, if `--sf=` points to a different file,
you need to make sure that it can be found through `INPUTPATH` in your configuration file.

## To conduct a heuristic algorithm

As mentioned in one of the arguments, the tool is able to run a number of built in heuristic algorithm for making adaptations.
The heuristics will also uses MIP solvers but with a much chaper computational costs. Currently, this repo provide the following heuristics:

* reactor
* highland
* bathtub
* extreme

It is flexible to run these heuristics with a specific statistical metric on the scenario inputs. Here is an example,
```julia
> julia adcc.jl ieee118 dc heuristic --sm=file-90-perc --sf="./format/scenario_example.json" --T=5 --EPS=0.0 --pf=paper_param --SOLVER=Cplex --HEURISTIC=reactor
```
The above command runs the algorithm using `reactor` heuristic on the 90-th percentile scenario of scenario file input `scenario_example.json`.

## To evaluate a design
Usually, we wish to see how a certain desig n performs in the future stochasticity. However, since the future is anticipative, we simulate, given this design, "how bad" (i.e. feasibility/load shedding over time) that future would be. To do so, you will need to store the design in an .json file just like the example in `./format/design.json`. Then, all you need to do is,
```julia
> julia adcc.jl ieee14 network evaluate --df="./format/design.json" --EVALOBJ=feasibility --S=100 --STOCHMODE=evolving
```
This command will drive an evaluation on the indicated design using internally generated scenarios under mode "evolving". The evaluation will take feasibility as the objective metric (--EVALOBJ).
If you wish to run the same evaluation with external scenarios,
```julia
> julia adcc.jl ieee14 network evaluate --DESIGNFILE="./format/design.json" --EVALOBJ=feasibility --STOCHFILE="./format/stoc.json"
```


## All command line arguments
The optimization model provides a variety of command line inputs for analyses. Here is a full list of them,
* Positional Arguments
    * `PROBLEM` (default=`ieee118`) problem network
        * link to `.m` file in `instance` dir
    * `MODEL` (default=`dc`) dispatch mode
        * `capacity` | `network`
    * `ALGORITHM` (default=`regular`) algorithm used for solving optimization
        * `sbd_norisk` | `sbd_heuristic`
* Keyword Arguments
    * `--T=1` (int) model total time periods, must comply with parameters and scenario inputs
    * `--EPS=0.0` (float) risk parameter
    * `--sm=file` (str) stochastic mode, give `file` for `.json` inputs
        * Can use post-fix terms to require specific statistical estimat
    * `--sf=` (str) is `--sm=file`, then must specify the `.json` path
    * `--df=` (str) points to the solution `.json` file (or dir) that needs to be evaluated for quality
    * `--eo=feasibility` (str) metric used the evaluation
    * `--COSTLambda=-1.0` (float) [0,1] parameter to skew the cost balance between expansion and hardening
    * `--DEMANDLambda=0.03` (float) percentage of demand growth in each time period
    * `--SHEDLambda=0.95` (float) minimum percentage of the load must be satisfied for feasibility Criteria
    * `--CONGESTLambda=1.0` (float) ratio of the original network congestion (thermal limits) to apply
    * `--ANGLESHIFTLambda=30` (float) degree of voltage angle allowed between two linked buses
    * `--DISCOUNTLambda=-0.01` (float) degragation rate of cost of adptation over time periods
    * `--SOLVER=config` (str) MIP solver selection, can also be `Cplex` or `Gurobi`. `config` means comply with config file, same as all the followings
    * `--PARALLEL=config` (str) whether or not to run `sbd_heuristic` or `sbd_norisk` in parallel.
    * `--WORKERS=config` (str) how many workers is required for parallel algorithm
    * `--TIMELIMIT=-1` (float) total time limit on running the optimization, `-1` means comply with config file, same as all the followings
    * `--TIMELIMITII=-1` (float) total time for finer control of the SBD algorithm
    * `--CGHEURISTIC=improver_heu` (str) heuristic algorithm used in SBD algorithm
    * `--USESBDNORIK=config` (bool) finer SBD algorithm control for better performance
    * `--HEURISTIC=reactor` (str) heuristic algorithm for optmization
    * `--NAME=00` (str) ouput name
    * `--OUTPUTPATH=` (str) additional output path controller for large scale testing

## Help
Feel free to explore using the following command,
```julia
> julia adcc.jl --help
```

## Problem Instances
Please find the power in the instance folder.

[1]:https://github.com/JuliaOpt/Gurobi.jl
[2]:https://github.com/JuliaOpt/CPLEX.jl
