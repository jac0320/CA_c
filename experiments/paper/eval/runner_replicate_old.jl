function write_batch_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2650_v2\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
end

function write_batch_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#PBS -N $(jobname)\n")
    write(f, "#PBS -m abe\n")
    wirte(f, "#PBS -M siteresults@gmail.com\n")
    write(f, "#PBS -l select=1:ncpus=2:mem=32gb,walltime=00:20:00\n")
    write(f, "#PBS -o /home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname)_o.out\n")
    write(f, "#PBS -e /home/sitew/Outputs/CLIMATE-PAPER/EVAL/$(outname)_e.out\n\n\n")
    write(f, "")
end


modegroup = ["network", "dc"]
subdirs = ["REACTOR/", "BATHTUB/", "HIGHLAND/"]

for sub in subdirs
    for mode in modegroup
        for seed in 55:60
        	filename = "paper_eval_replicate_s$(seed)_1000_$(mode).sh"
            inputname = "paper_s$(seed)_1000"
        	jobname = "p_eval_s$(seed)_1000_$(mode)"
        	outname = "paper_eval_replicate_s$(seed)_1000_$(mode)"
        	if !isfile(filename)
            	f = open(filename, "w")
            	write_batch_head(f, jobname, outname)
        		write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 $(mode) evaluate --TIMELIMIT=86400 --pf=paper_param --sm=file --sf=$(inputname).json --eo=feasibility --df=/home/sitew/Outputs/CLIMATE-PAPER/REPLICATE_DESIGN/$(sub) --T=5 --SOLVER=Gurobi")
            	close(f)
            	run(`sbatch $(filename)`)
        	end
        end
    end
end
