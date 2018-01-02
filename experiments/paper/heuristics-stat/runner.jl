function write_slurm_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/HEURISTIC-STAT/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2660_v3\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
    # write(f, "julia /home/sitew/start_gurobi.jl\n")
end

function write_pbs_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#PBS -N $(jobname)\n")
    write(f, "#PBS -l select=1:ncpus=2:mem=4gb,walltime=00:20:00\n")
    write(f, "#PBS -o /home/sitew/Outputs/CLIMATE-PAPER/HEURISTIC-STAT/$(outname)_o.out\n")
    write(f, "#PBS -e /home/sitew/Outputs/CLIMATE-PAPER/HEURISTIC-STAT/$(outname)_e.out\n\n\n")
end

sizegroup = [50, 100, 200, 500, 1000, 2000]
modegroup = ["dc"]
heuristicgroup = ["reactor", "highland", "bathtub"]
statgroup = ["-90-perc", "-95-perc", "-max", "-average"]

seedtotal_s = 1
seedtotal_e = 10

for mode in modegroup
    for size in sizegroup
        for seed in seedtotal_s:seedtotal_e
            for heu in heuristicgroup
                for stat in statgroup
            		filename = "p_heustat_$(heu)$(stat)_s$(seed)_$(size)_$(mode).sh"
                    inputname = "paper_s$(seed)_$(size)"
            		jobname = "p_heustat_$(heu)$(stat)_s$(seed)_$(size)_0_$(mode)"
            		outname = jobname
                    f = open(filename, "w")
                    if ARGS[1] == "slurm"
                        write_slurm_head(f, jobname, outname)
                        write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) heuristic --pf=paper_param --sm=file$(stat) --sf=$(inputname).json --HEURISTIC=$(heu) --NAME=$(outname) --T=5 --OUTPUTPATH=HEUSTAT_DESIGN/ --SOLVER=CPLEX")
                    elseif ARGS[1] == "pbs"
                        write_pbs_head(f, jobname, outname)
                        write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) heuristic --pf=paper_param --sm=file$(stat) --sf=$(inputname).json --HEURISTIC=$(heu) --NAME=$(outname) --T=5 --OUTPUTPATH=HEUSTAT_DESIGN/ --SOLVER=Gurobi")
                    end
                    close(f)
                    run(`sbatch $(filename)`)
                    rm(filename)
                end
            end
        end
    end
end
