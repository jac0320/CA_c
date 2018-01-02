#!/bin/bash
#SBATCH --job-name="rcp85_100"
#SBATCH --output="/home/sitew/Outputs/Climate-PAPER/WORKERS_SBD8_100_20.out"
#SBATCH --constraint="cpu_model:E5-2699_v4"#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/Github/climate/adcc.jl ieee118 dc sbd_heuristic --PARALLEL=1 --WORKERS=8 --CGHEURISTIC=lineup_partial_heuristic --SOLVER=CPLEX --sm=file --sf=rcp85_100.json --EPS=0.2 --NAME=WORKERS_SBD8_100_20