function write_slurm_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/EXP4/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2660_v3\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 144:00:00\n\n")
end

function write_pbs_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#PBS -N $(jobname)\n")
    write(f, "#PBS -m abe\n")
    write(f, "#PBS -M siteresults@gmail.com\n")
    write(f, "#PBS -l select=1:ncpus=16:mem=64gb,walltime=72:00:00\n")
    write(f, "#PBS -o /home/sitew/Outputs/CLIMATE-PAPER/RISK/$(outname).o\n")
    write(f, "#PBS -e /home/sitew/Outputs/CLIMATE-PAPER/RISK/$(outname).e\n\n\n")
    write(f, "module add gurobi/7.0.1 \n")
end

slrgroup = [1, 25, 50, 75, 98]
slrgroup = [5,10,15,20,30,35,40,45,55,60,65,70,80,85,90,95]
sizegroup = [100, 200]
epsgroup = [0.0, 0.05]
epsgroup = [0.0]
start_seed = 3
end_seed = 5

algmap = Dict(true=>"sbd_norisk", false=>"sbd_heuristic")

for size in sizegroup
    for slr in slrgroup
    	for seed in start_seed:end_seed
            for eps in epsgroup
                inputname = "paper2_slr$(slr)_s$(seed)_$(size)"
        		jobname = "p2_exp4_slr$(slr)_s$(seed)_$(size)_$(Int(eps*100))_dc"
        		outname = "p2_exp4_slr$(slr)_s$(seed)_$(size)_$(Int(eps*100))_dc"
                if ARGS[1] == "slurm"
                    filename = "p2_exp4_slr$(slr)_s$(seed)_$(size)_$(Int(eps*100))_dc.sh"
                    f = open(filename, "w")
                    write_slurm_head(f, jobname, outname)
                elseif ARGS[1] == "pbs"
                    filename = "p2_exp4_slr$(slr)_s$(seed)_$(size)_$(Int(eps*100))_dc.pbs"
                    f = open(filename, "w")
                    write_pbs_head(f, jobname, outname)
                end
    			write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 dc $(algmap[eps==0.0]) --TIMELIMITII=3600 --SOLVER=CPLEX --CGHEURISTIC=lineup_partial_heuristic --pf=paper_param --sm=file --USESBDNORISK=1 --PARALLEL=1 --WORKERS=10 --sf=$(inputname).json --EPS=$(eps) --NAME=$(outname) --T=5")
                close(f)
                if ARGS[1] == "slurm"
                    run(`sbatch $(filename)`)
                elseif ARGS[1] == "pbs"
                    run(`qsub $(filename)`)
                end
                rm(filename)
            end
        end
    end
end
