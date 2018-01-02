#!/bin/bash
#SBATCH --job-name="rcp85_20"
#SBATCH --output="/home/sitew/Outputs/Climate/BENCH_85_20_0.out"
#SBATCH --constraint="cpu_model:E5-2660_v3"
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/Github/sbd-climate/adcc.jl ieee118 dc regular --SOLVER=CPLEX --STOCHMODE=file --TIMELIMIT=21600 --STOCHFILE=rcp85_20.json --EPS=0.0 --NAME=BENCH_85_20_0