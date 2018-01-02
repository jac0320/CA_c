sizegroup = [200, 300, 500, 1000, 2000]
modegroup = ["dc"]
heuristicgroup = ["reactor", "highland", "bathtub"]
statgroup = ["-90-perc", "-95-perc", "-max", "-average"]

seedtotal_s = 1
seedtotal_e = 5

f = open("heustat_local.sh", "w")

for mode in modegroup
    for size in sizegroup
        for seed in seedtotal_s:seedtotal_e
            for heu in heuristicgroup
                for stat in statgroup
                    inputname = "paper_s$(seed)_$(size)"
            		outname = "p_heustat_$(heu)$(stat)_s$(seed)_$(size)_0_$(mode)"
                    if is_windows()
                        write(f, "julia C:\\\\Users\\\\bitja\\\\Dropbox\\\\SourceCodes\\\\ADCC\\\\adcc.jl ieee118 $(mode) heuristic --pf=paper_param --sm=file$(stat) --sf=$(inputname).json --HEURISTIC=$(heu) --NAME=$(outname) --T=5 --OUTPUTPATH=HEUSTAT_DESIGN/ --SOLVER=Gurobi\n")
                    elseif is_apple() || is_linux()
                        write(f, "julia $(homedir())/Dropbox/SourceCodes/ADCC/adcc.jl ieee118 $(mode) heuristic --pf=paper_param --sm=file$(stat) --sf=$(inputname).json --HEURISTIC=$(heu) --NAME=$(outname) --T=5 --OUTPUTPATH=HEUSTAT_DESIGN/ --SOLVER=Cplex\n")
                    end
                end
            end
        end
    end
end
