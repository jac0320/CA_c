#!/bin/bash
#SBATCH --job-name="rcp85_100"
#SBATCH --output="/home/sitew/Outputs/Climate/BENCH_85_100_80.out"
#SBATCH --constraint="cpu_model:E5-2660_v3"
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/Github/sbd-climate/adcc.jl ieee118 dc regular --SOLVER=CPLEX --STOCHMODE=file --TIMELIMIT=21600 --STOCHFILE=rcp85_100.json --EPS=0.8 --NAME=BENCH_85_100_80