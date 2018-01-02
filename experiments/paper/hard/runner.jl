function write_batch_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/HARD/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2695_v4\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 25:00:00\n\n\n")
    write(f, "julia /home/sitew/start_gurobi.jl\n\n")
end


function write_pbs_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#PBS -N $(jobname)\n")
    write(f, "#PBS -m abe\n")
    write(f, "#PBS -M siteresults@gmail.com\n")
    write(f, "#PBS -l select=1:ncpus=16:mem=64gb,walltime=72:00:00\n")
    write(f, "#PBS -o /home/sitew/Outputs/CLIMATE-PAPER/HARD/$(outname)_o.out\n")
    write(f, "#PBS -e /home/sitew/Outputs/CLIMATE-PAPER/HARD/$(outname)_e.out\n\n\n")
    write(f, "module add gurobi/7.0.1 \n")
end

sizegroup = [20, 50, 100, 200]
modegroup = ["network", "dc"]
epsgroup = [0.0, 0.05, 0.1, 0.2, 0.3]
seedtotal = 5

for size in sizegroup
    for mode in modegroup
        for seed in 1:seedtotal
            for eps in epsgroup
                inputname = "paper_s$(seed)_$(size)"
        		jobname = "p_hard_s$(seed)_$(size)_$(Int(eps*100))_$(mode)"
        		outname = jobname
                if ARGS[1] == "slurm"
                    filename = "p_hard_s$(seed)_$(size)_$(Int(eps*100))_$(mode).sh"
                    f = open(filename, "w")
                    write_slurm_head(f, jobname, outname)
                else
                    filename = "p_hard_s$(seed)_$(size)_$(Int(eps*100))_$(mode).pbs"
                    f = open(filename, "w")
                    write_pbs_head(f, jobname, outname)
                end
    			write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) regular --TIMELIMIT=172800 --SOLVER=Gurobi --pf=paper_param --sm=file --sf=$(inputname).json --EPS=$(eps) --NAME=$(outname) --OUTPUTPATH=HARD_DESIGN/ --T=5")
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
