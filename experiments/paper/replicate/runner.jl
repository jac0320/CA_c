function write_slurm_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/REPLICATE-v2/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2695_v4\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
    write(f, "julia /home/sitew/start_gurobi.jl\n")
end

function write_pbs_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#PBS -N $(jobname)\n")
    write(f, "#PBS -m abe\n")
    write(f, "#PBS -M siteresults@gmail.com\n")
    write(f, "#PBS -l select=1:ncpus=16:mem=64gb,walltime=72:00:00\n")
    write(f, "#PBS -o /home/sitew/Outputs/CLIMATE-PAPER/REPLICATE/$(outname).o\n")
    write(f, "#PBS -e /home/sitew/Outputs/CLIMATE-PAPER/REPLICATE/$(outname).e\n\n\n")
    write(f, "module add gurobi/7.0.1 \n")
end

sizegroup = [200, 300, 400, 500]
modegroup = ["capacity", "network", "dc"]
seedtotal = 10
algmap = Dict("network"=>"sbd_norisk", "dc"=>"sbd_norisk", "capacity"=>"regular")

for mode in modegroup
    for size in sizegroup
        for seed in 1:seedtotal
            inputname = "paper2_s$(seed)_$(size)"
    		jobname = "p2_rnr_s$(seed)_$(size)_0_$(mode)"
    		outname = "p2_rnr_s$(seed)_$(size)_0_$(mode)"
            if ARGS[1] == "slurm"
                filename = "p2_rnr_s$(seed)_$(size)_0_$(mode).sh"
                f = open(filename, "w")
                write_slurm_head(f, jobname, outname)
            elseif ARGS[1] == "pbs"
                filename = "p2_rnr_s$(seed)_$(size)_0_$(mode).pbs"
                f = open(filename, "w")
                write_pbs_head(f, jobname, outname)
            end
			write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) $(algmap[mode]) --SOLVER=Gurobi --TIMELIMITII=7200 --pf=paper_param --sm=file --sf=$(inputname).json --EPS=0.0 --NAME=$(outname) --OUTPUTPATH=REPLICATE_DESIGN/ --T=5")
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
