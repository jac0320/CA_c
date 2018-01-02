#!/bin/bash
#SBATCH --job-name="rcp85_150"
#SBATCH --output="/home/sitew/Outputs/Climate/SBD_85_150_70_PA.out"
#SBATCH --constraint="cpu_model:E5-2660_v3"#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/Github/sbd-climate/adcc.jl ieee118 dc sbd_heuristic --PARALLEL=1 --WORKERS=16 --SOLVER=CPLEX --STOCHMODE=file --CGHEURISTIC=lineup_partial_heuristic --STOCHFILE=rcp85_150.json --EPS=0.7 --NAME=SBD_85_150_70_PA