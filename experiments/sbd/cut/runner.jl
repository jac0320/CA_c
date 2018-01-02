function write_batch_head(f, jobname, outname)
    write(f, "#!/bin/bash\n")
    write(f, "#SBATCH --job-name=\"$(jobname)\"\n")
    write(f, "#SBATCH --output=\"/home/sitew/Outputs/Climate/$(outname).out\"\n")
    write(f, "#SBATCH --constraint=\"cpu_model:E5-2660_v3\"")
    write(f, "#SBATCH --ntasks=1\n")
    write(f, "#SBATCH --no-requeue\n")
    write(f, "#SBATCH -t 72:00:00\n\n")
    # write(f, "julia /home/sitew/start_gurobi.jl\n")
end

probgroup = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
sizegroup = [20, 30, 40, 50, 60, 80]

cggroup = ["", "--CGHEURISTIC=lineup_partial_heuristic"]
cgname = Dict(""=>"IM", "--CGHEURISTIC=lineup_partial_heuristic"=>"PA")

for s in sizegroup
    for p in probgroup
        for cg in cggroup
            fname = "rcp85_$(s)_$(Int(p*100))_$(cgname[cg]).sh"
            pname = "--EPS=$(p)"
            jname = "rcp85_$(s)"
            oname = "FIX_SBD_85_$(s)_$(Int(p*100))_$(cgname[cg])"
            f = open(fname, "w")
            write_batch_head(f, jname, oname)
            write(f, "julia /home/sitew/Github/climate/adcc.jl ieee118 dc sbd_heuristic --PARALLEL=1 --WORKERS=16 --SOLVER=CPLEX --sm=file $(cg) --sf=$(jname).json $(pname) --NAME=$(oname)")
            close(f)
            run(`sbatch $(fname)`)
        end
    end
end
