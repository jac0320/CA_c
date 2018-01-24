function write_output_files(power::Dict, param::Dict, stoc::stocType, solution, driver::Dict)

	designFilepath = string(driver[:OUTPUTPATH], "design_", driver[:NAME],".json")
	paramFilepath = string(driver[:OUTPUTPATH], "param_", driver[:NAME],".json")
	stocFilepath = string(driver[:OUTPUTPATH], "stoc_",stoc.S,"_",driver[:STOCHMODE],"_",driver[:NAME],".json")

	write_json(power,designFilepath)

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

	return
end

function summary_driver_arguments(param::Dict, stoc::stocType, driver::Dict)

    info("Problem Instance      : ", driver[:PROBLEM])
    info("Characteristic        : ", driver[:MODEL])
    info("Stochastic Mode       : ", driver[:STOCHMODE])
    info("Algorithm             : ", driver[:ALGO])
    info("Time Periods(T)       : ", driver[:T])
    info("Scenario Count(S)     : ", driver[:S])
    info("Risk(eps)             : ", driver[:eps])
	info("Demand Change         : ", driver[:DEMANDLambda])
    info("Shedding Allowing     : ", driver[:SHEDLambda])
	info("Congestion            : ", driver[:CONGESTLambda])
	info("Angle Shift Limit     : ", driver[:ANGLESHIFTLambda])
	info("Discounting Cost      : ", driver[:DISCOUNTLambda])
	info("Cost Ratio (expand)   : ", driver[:COSTLambda])
	info("Time Limit L1         : ", driver[:TIMELIMIT])
	info("Time Limit L2         : ", driver[:TIMELIMITII])
	info("Time Limit L3         : ", driver[:TIMELIMITIII])
	info("Parallel Indicator    : ", driver[:PARALLEL])
	info("Workers Utilized      : ", driver[:WORKERS])
    info("Single Worker Threads : ", driver[:WORKERTHREADS])
	info("Warm Start            : ", driver[:WARMSTART])
	info("Job created at        : ", now())
	info("Job Output Name       : ", driver[:NAME])
	info("Julia Version         : ", VERSION)
	info("CPU Cores             : ", Sys.CPU_CORES)
	info("Machine INFO          : ", Sys.MACHINE)
	info("CPU Summary           : ")
	Sys.cpu_summary()

	return
end

function write_json(content::Dict,filename::AbstractString, pathprefix::AbstractString="")

	isempty(pathprefix) ? wf = open(filename, "w") : wf = open(pathprefix,"/",filename,"w")

	write(wf, JSON.json(content))
	close(wf)

	return
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

function get_cache_code(algo::AbstractString, stage::AbstractString, driver::Dict)

	cacheCode = "./experiments/cache/"

	cacheCode = string(cacheCode, driver[:PROBLEM], driver[:MODEL], replace(driver[:STOCHFILE],".json",""))
	cacheCode = string(cacheCode, algo, stage)

	return cacheCode

end

function slr_shutdown(param::Dict, scenario::scenarioType)

	B, T = param[:B], param[:T]

	cP = ones(B, T)
	for b in 1:B
		for t in 1:T
			scenario.data["SL"][t] >= param[:Ele][b] ? cP[b,t] = 0 : cP[b,t] = 1
		end
	end

	return cP
end

function static_expansion(param::Dict, driver::Dict; hl=nothing)

	T, B = param[:T], param[:B]

	zeroStoc = null_scenario_stoc(B, T)
	eP = build_sp(param, zeroStoc, driver, selection=[1], sbtype="tight")

	if hl != nothing
		for i=1:B
			param[:Ele][i] < hl && @constraint(eP.model, [t=1:T], eP.vars[:pg][i,t] <= param[:Pg0][i])
		end
	end

	config_solver(eP.model, driver, mipgap=0.0001)
	status = solve(eP.model, suppress_warnings=true)

	status == :Infeasible && print_iis_gurobi(eP.model, driver)

	return getvalue(eP.vars[:pg])

end

function threshold_harden(param::Dict, scen::scenarioType, threshold::Float64, lag::Bool, closeP::Any)

	B, T = param[:B], param[:T]

	hP = zeros(B, T)

	if lag
		for t in 2:T, b in 1:B
			hP[b,t] = hP[b,t-1]
			if closeP[b,t] == 1 && scen.data["SS"][b,t-1] > param[:Ele][b]
				hP[b,t] = max(hP[b,t-1], ceil((scen.data["SS"][b,t-1] - param[:Ele][b])/param[:ProM][b]))
			end
		end
		for b in 1:B
			hP[b, 1] = floor(hP[b, 2]/2)
		end
	else
		for t in 1:T, b in 1:B
			if t >= 2
				hP[b,t] = hP[b,t-1]
			end
			if closeP[b,t] == 1 && scen.data["SS"][b,t] >= param[:Ele][b]
				if t > 1
					hP[b,t] = max(hP[b,t-1], ceil((threshold - param[:Ele][b])/param[:ProM][b]))
				else
					hP[b,t] = ceil((threshold - param[:Ele][b])/param[:ProM][b])
				end
			end
		end
	end

	return hP
end

function enforce_bound(prob::oneProblem, varkey::Symbol;lb::Vector=[], ub::Vector=[])

	if !isempty(lb)
		@assert size(prob.vars[varkey]) == size(lb)
		for (i,k) in enumerate(prob.vars[varkey])
			setlowerbound(k, lb[i])
		end
	end

	if !isempty(ub)
		@assert size(prob.vars[varkey]) == size(ub)
		for (i,k) in prob.vars[varkey]
			setupperbound(k, ub[i])
		end
	end

	return
end

function enforce_startval(prob::oneProblem, varkey::Symbol; sval::Vector=[])

	isempty(sval) && return

	@assert size(sval) == size(prob.vars[varkey])

	for (i,k) in prob.vars[varkey]
		setvalue(k, sval[i])
	end

	return
end

function enforce_startval(prob::oneProblem; sval::Vector=[])
	isempty(sval) && return

	@assert length(sval) == length(prob.model.colVal)

	C = length(prob.model.colVal)
	for i in 1:C
		setvalue(Variable(prob.model, i), sval[i])
	end

	return
end
