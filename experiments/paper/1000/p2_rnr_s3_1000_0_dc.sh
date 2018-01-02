#!/bin/bash
#SBATCH --job-name="p2_rnr_s3_1000_0_dc"
#SBATCH --output="/home/sitew/Outputs/CLIMATE-PAPER/REPLICATE/p2_rnr_s3_1000_0_dc.out"
#SBATCH --constraint="cpu_model:E5-2660_v3"
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/start_gurobi.jl
julia /home/sitew/Github/climate/adcc.jl ieee118 dc sbd_norisk --SOLVER=Gurobi --TIMELIMITII=3600 --pf=paper_param --sm=file --sf=paper2_s3_1000.json --EPS=0.0 --NAME=p2_rnr_s3_1000_0_dc --OUTPUTPATH=REPLICATE_DESIGN/ --T=5