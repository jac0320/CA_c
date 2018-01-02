function write_slurm_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2650_v2\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
end

function write_pbs_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#PBS -N $(jobname)\n")
    write(f, "#PBS -m abe\n")
    write(f, "#PBS -M siteresults@gmail.com\n")
    write(f, "#PBS -l select=1:ncpus=8:mem=32gb,walltime=72:00:00\n")
    write(f, "#PBS -o /home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname)_o.out\n")
    write(f, "#PBS -e /home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname)_e.out\n\n\n")
end


modegroup = ["network", "dc"]
subdirs = ["REACTOR/", "BATHTUB/", "HIGHLAND/"]

for sub in subdirs
    for mode in modegroup
        for seed in 55:60
            inputname = "paper_s$(seed)_1000"
        	jobname = "p_eval-$(sub[1:end-1])_s$(seed)_1000_$(mode)"
        	outname = "p_eval-$(sub[1:end-1])_s$(seed)_1000_$(mode)"
            if ARGS[1] == "slurm"
                filename = "p_eval-$(sub[1:end-1])_s$(seed)_1000_$(mode).sh"
        	    f = open(filename, "w")
        	    write_batch_head(f, jobname, outname)
            elseif ARGS[1] == "pbs"
                filename = "p_eval-$(sub[1:end-1])_s$(seed)_1000_$(mode).pbs"
                f = open(filename, "w")
                write_pbs_head(f, jobname, outname)
            end
    		write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) evaluate --pf=paper_param --sm=file --sf=$(inputname).json --eo=feasibility --df=/home/sitew/Outputs/CLIMATE-PAPER/HEURISTIC-STAT/$(sub) --NAME=$(outname) --T=5 --SOLVER=Cplex")
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
