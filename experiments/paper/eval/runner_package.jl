function write_slurm_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname).out\"\n")
    # write(f, "#SBATCH --constraint=\"cpu_model:E5-2650_v2\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
end

function write_pbs_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#PBS -N $(jobname)\n")
    write(f, "#PBS -m abe\n")
    write(f, "#PBS -M sitew@g.clemson.edu\n")
    write(f, "#PBS -l select=1:ncpus=8:mem=32gb,walltime=72:00:00\n")
    write(f, "#PBS -o /home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname)_o.out\n")
    write(f, "#PBS -e /home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname)_e.out\n\n\n")
end

modegroup = ["network", "dc"]
sizegroup = [500, 1000]
sizegroup = [2000, 5000]

for size in sizegroup
    for mode in modegroup
        for seed in 1:10
            inputname = "paper_s$(seed)_$(size)"
        	jobname = "p_package_eval_s$(seed)_$(size)_$(mode)"
        	outname = "p_package_eval_s$(seed)_$(size)_$(mode)"
            if ARGS[1] == "slurm"
                filename = "p_package_eval_s$(seed)_$(size)_$(mode).sh"
        	    f = open(filename, "w")
        	    write_slurm_head(f, jobname, outname)
            elseif ARGS[1] == "pbs"
                filename = "p_package_eval_s$(seed)_$(size)_$(mode).pbs"
                f = open(filename, "w")
                write_pbs_head(f, jobname, outname)
            end
    		write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) evaluate --pf=paper_param --sm=file --sf=$(inputname).json --eo=feasibility --df=/home/sitew/Outputs/CLIMATE-PAPER/EVAL_PACKAGE/ --NAME=$(outname) --T=5 --SOLVER=Cplex")
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
