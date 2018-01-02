#!/bin/bash
#SBATCH --job-name="p_risk_s3_1000_10_dc"
#SBATCH --output="/home/sitew/Outputs/CLIMATE-PAPER/RISK/p_risk_s3_1000_10_dc.out"
#SBATCH --constraint="cpu_model:E5-2660_v3"
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH -t 72:00:00

julia /home/sitew/Github/climate/adcc.jl ieee118 dc sbd_heuristic --TIMELIMITII=3600 --SOLVER=CPLEX --CGHEURISTIC=lineup_partial_heuristic --pf=paper_param --sm=file --USESBDNORISK=1 --PARALLEL=1 --WORKERS=4 --sf=paper_s3_1000.json --EPS=0.1 --NAME=p_risk_s3_1000_10_dc --T=5 --OUTPUTPATH=RISK/