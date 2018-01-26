# Energy Infrastructure Adaptive Design for Climate Change

This project is designed for exploring resilient power systems designs in response to future climate changes, which is my dissertation all about.

## To do a adaptive design optimization
```julia
> julia adcc.jl ieee118 dc regular --sm=file --sf="./format/scenario_example.json" --T=5 --EPS=0.0 --pf=paper_param --SOLVER=Cplex --TIMELIMIT=300 --NAME=yourscenario
```

## To conduct a heuristic algorithm
```julia
> julia adcc.jl ieee118 dc heuristic --sm=file-90-perc --sf="./format/scenario_example.json" --T=5 --EPS=0.0 --pf=paper_param --SOLVER=Cplex --HEURISTIC=reactor
```
The above command runs the algorithm using `reactor` heuristic on the 90-th percentile scenario of scenario file input `scenario_example.json`.

## To evaluate a design
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
        * `shcg_nr` | `shcg`
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
    * `--PARALLEL=config` (str) whether or not to run `shcg` or `shcg_nr` in parallel.
    * `--WORKERS=config` (str) how many workers is required for parallel algorithm
    * `--TIMELIMIT=-1` (float) total time limit on running the optimization, `-1` means comply with config file, same as all the followings
    * `--TIMELIMITII=-1` (float) total time for finer control of the SBD algorithm
    * `--CGHEURISTIC=improver_heu` (str) heuristic algorithm used in SBD algorithm
    * `--USESBDNORIK=config` (bool) finer SBD algorithm control for better performance
    * `--HEURISTIC=reactor` (str) heuristic algorithm for optmization
    * `--NAME=00` (str) ouput name
    * `--OUTPUTPATH=` (str) additional output path controller for large scale testing
