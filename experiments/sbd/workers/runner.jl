function write_batch_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/Climate-PAPER/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2699_v4\"")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
    # write(f, "julia /home/sitew/start_gurobi.jl\n")
end

probgroup = [0.05, 0.1, 0.2, 0.3]
workergroup = [2,4,8]


    for p in probgroup
      	for w in workergroup
            fname = "rcp85_100_$(Int(p*100)).sh"
            pname = "--EPS=$(p)"
            jname = "rcp85_100"
			oname = "WORKERS_SBD$(w)_100_$(Int(p*100))"
            f = open(fname, "w")
            write_batch_head(f, jname, oname)
			write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 dc sbd_heuristic --PARALLEL=1 --WORKERS=$(w) --CGHEURISTIC=lineup_partial_heuristic --SOLVER=CPLEX --sm=file --sf=rcp85_100.json $(pname) --NAME=$(oname)")
            close(f)
            run(`sbatch $(fname)`)
        end
    end
