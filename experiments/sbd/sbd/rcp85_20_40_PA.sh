#!/bin/bash
#SBATCH --job-name="rcp85_20"
#SBATCH --output="/home/sitew/Outputs/Climate/FIX_SBD_85_20_40_PA.out"
#SBATCH --constraint="cpu_model:E5-2660_v3"#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/Github/climate/adcc.jl ieee118 dc sbd_heuristic --PARALLEL=1 --WORKERS=16 --SOLVER=CPLEX --sm=file --CGHEURISTIC=lineup_partial_heuristic --sf=rcp85_20.json --EPS=0.4 --NAME=FIX_SBD_85_20_40_PA