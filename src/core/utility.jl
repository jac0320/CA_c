function write_output_files(power::Dict, param::Dict, stoc::stocType, solution, exargs::Dict)

	if config.FILEOUTPUT == true

		# Assign File Output Strings
		designFilepath = string(config.OUTPUTPATH, exargs[:OUTPATH], "design_", exargs[:NAME],".json")
		paramFilepath = string(config.OUTPUTPATH, exargs[:OUTPATH], "param_", exargs[:NAME],".json")
		stocFilepath = string(config.OUTPUTPATH, exargs[:OUTPATH], "stoc_",stoc.S,"_",exargs[:STOCHMODE],"_",exargs[:NAME],".json")

		# See if there already exist such a file ::
		try
			rm(designFilepath)
		catch e
			info("Creating new design output file.")
		end
		write_json(power,designFilepath)

		if exargs[:PARAMWRITE]
			try
				rm(paramFilepath)
			catch e
				info("Creating new param output file")
			end
			write_json(param,paramFilepath)
		end

		# Write scnearios down
		if exargs[:STOCHWRITE]
			try
				rm(stocFilePath)
			catch e
				println("Creating new stoc output file")
			end
			stocDict = Dict()
			stocDict["SS"] = Dict()
			stocDict["SL"] = Dict()
			for s in 1:stoc.S
				stocDict["SS"] = stoc.scenarios[s].data["SS"]
				stocDict["SL"] = stoc.scenarios[s].data["SL"]
			end
			write_json(stocDict,stocFilepath)
		end

		# Write Design Details
		designDict = Dict()
        if isa(solution, solnType)
			designDict["pg"] = solution.primal[:pg]
			designDict["h"] = solution.primal[:h]
		elseif isa(solution, designType)
			designDict["pg"] = solution.pg
			designDict["h"] = solution.h
		else
			error("Unkown solution data structure.")
		end
		write_json(designDict, designFilepath)
	end
end

function screen_output(solution::solnType)
	if config.SCREENSHOW == true
		info(":: Solution pg :: ")
		println(solution.primal[:pg])
		info(":: Solution h :: ")
		println(solution.primal[:h])
		flush(STDOUT)
	end
end

function write_json(content::Dict,filename::AbstractString, pathprefix::AbstractString="")

	# File name must have .suffix
	if isempty(pathprefix)
		wf = open(filename, "w")
	else
		wf = open(pathprefix,"/",filename,"w")
	end

	write(wf, JSON.json(content))
	close(wf)

end


function arrarr2mat(arrarr::Array, numType, dimA::Int, dimB::Int)

	@assert dimB == length(arrarr)
	@assert dimA == length(arrarr[1])

	mat = zeros(numType, dimA, dimB)
	for b in 1:dimB
		for a in 1:dimA
			mat[a,b] = arrarr[a][b]
		end
	end

	return mat
end

function get_rhs(model::JuMP.Model)

    rowCnt = MathProgBase.numlinconstr(model)
    rhs = Array{Float64}(rowCnt)

    for (nrow, con) in enumerate(model.linconstr)
        # There shouldn't be any two sided constraints
        if (con.lb == 0 && con.ub == Inf) || (con.lb == -Inf && con.ub == 0) || (con.lb == 0 && con.ub == 0)
            rhs[nrow] = 0.0
        elseif (con.lb == -Inf) # >= constraint with non-zero rhs
            rhs[nrow] = con.ub
        elseif (con.ub == Inf)  # >= constraint with non-zero rhs
            rhs[nrow] = con.lb
        else                    # Equality constraint with non-zero rhs
            rhs[nrow] = conlb
        end
    end

    return rhs
end

function get_cache_code(algo::AbstractString, stage::AbstractString, exargs::Dict)

	cacheCode = "./experiments/cache/"

	cacheCode = string(cacheCode, exargs[:PROBLEM], exargs[:MODEL], replace(exargs[:STOCHFILE],".json",""))
	cacheCode = string(cacheCode, algo, stage)

	return cacheCode

end
