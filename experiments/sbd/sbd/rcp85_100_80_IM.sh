#!/bin/bash
#SBATCH --job-name="rcp85_100"
#SBATCH --output="/home/sitew/Outputs/Climate/SBD_85_100_80_IM.out"
#SBATCH --constraint="cpu_model:E5-2660_v3"#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/Github/sbd-climate/adcc.jl ieee118 dc sbd_heuristic --PARALLEL=1 --WORKERS=16 --SOLVER=CPLEX --STOCHMODE=file  --STOCHFILE=rcp85_100.json --EPS=0.8 --NAME=SBD_85_100_80_IM