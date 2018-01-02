function write_batch_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/Climate/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2660_v3\"\n")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
    # write(f, "julia /home/sitew/start_gurobi.jl\n")
end

probgroup = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
sizegroup = [20, 50, 100, 150]

for s in sizegroup
    for p in probgroup
            fname = "rcp85_$(s)_$(Int(p*100)).sh"
            pname = "--EPS=$(p)"
            jname = "rcp85_$(s)"
			oname = "BENCH_85_$(s)_$(Int(p*100))"
            f = open(fname, "w")
            write_batch_head(f, jname, oname)
            write(f, "julia /home/sitew/Github/sbd-climate/adcc.jl ieee118 dc regular --SOLVER=CPLEX --STOCHMODE=file --TIMELIMIT=21600 --STOCHFILE=$(jname).json $(pname) --NAME=$(oname)")
            close(f)
            run(`sbatch $(fname)`)
    end
end
