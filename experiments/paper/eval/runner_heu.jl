function write_batch_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2660_v3\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
end

modegroup = ["network", "dc"]

for mode in modegroup
    for seed in 55:60
    	filename = "paper_eval_heu_s$(seed)_5000_$(mode).sh"
        inputname = "paper_s$(seed)_5000"
    	jobname = "p_eval_s$(seed)_5000_$(mode)"
    	outname = "paper_eval_heu_s$(seed)_5000_$(mode)"
    	if !isfile(filename)
        	f = open(filename, "w")
        	write_batch_head(f, jobname, outname)
    		write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) evaluate --pf=paper_param --sm=file --sf=$(inputname).json --eo=feasibility --df=/home/sitew/Outputs/CLIMATE-PAPER/HEURISTIC_DESIGN/ --T=5 --SOLVER=CPLEX")
        	close(f)
        	run(`sbatch $(filename)`)
    	end
    end
end
