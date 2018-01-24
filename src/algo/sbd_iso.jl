function isolate_stage(power::Dict, param::Dict, stoc::stocType, exargs::Dict,
					   isoprob_formulation::Function, selection=[];
					   kwargs...)

	options = Dict(kwargs)

	isoTime = 0.0
	skipsafety = false

	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	isempty(selection) ? S = stoc.S : S = length(selection)
	(isempty(selection)) && (selection = [1:S;])

	haskey(options, :skipsafety) ? skipsafety = true : skipsafety = false

	# =========================== Isolation Stage ===========================#
	isoDesignPool = Array{designType}(exargs[:S])
	isoDesignCost = Array{Float64}(exargs[:S])
	isoDesignCover = Array{Float64}(exargs[:S])
	isoDesignLB = Array{Float64}(exargs[:S])

	earlyExitIndex = 0
	earlyExitCost = Inf
	earlyExitFeaCnt = 0
	earlyExit = false

	if config.PARALLEL # parallel implementation
		info("[ISO] Running in parallel...")
		solveTime = @elapsed isoDesignPool = pmap((a1,a2,a3,a4,a5,a6,a7)->isolate_solve_one_scenario(a1,a2,a3,a4,a5,a6,a7),
							[power for s in selection],
							[param for s in selection],
							[stoc for s in selection],
							[s for s in selection],
							[exargs for s in selection],
							[isoprob_formulation for s in selection],
							[subprobType for s in selection])

		checkFeaTime = @elapsed isoDesignPool = pmap((a1,a2,a3,a4,a5,a6)->isolate_check_one_design_feasibility(a1,a2,a3,a4,a5,a6),
							[power for d in isoDesignPool],
							[param for s in isoDesignPool],
							[stoc for s in isoDesignPool],
							[exargs for d in isoDesignPool],
							[d for d in isoDesignPool],
							[S for d in isoDesignPool])

		for s in selection
			info("[ISO] Scen $s -> cost = [$(round.(isoDesignPool[s].cost,2))][LB=$(round.(isoDesignPool[s].lb,2))][TIME=$(round.(isoDesignPool[s].time,2))][Cover=$(round.(isoDesignPool[s].coverage,4))]")
			push!(stoc.sbdColumns, isoDesignPool[s])
		end

	else # sequential implementation
		info("[ISO] Running in sequential...")
		solveTime = 0.0
		for s in selection
			oneSolveTime = @elapsed isoDesign =
				isolate_solve_one_scenario(power, param, stoc, s, exargs, isoprob_formulation, subprobType)
			push!(stoc.sbdColumns, isoDesign)
			isoDesignPool[s] = isoDesign
			solveTime += oneSolveTime
		end

		info("[ISO] Enetering feasibility checking phase...")
		checkFeaTime = 0.0
		allSubprobs = Array{oneProblem}(length(selection))
		info("[ISO] WARMSTARTing feasibility check process...")
		for s in selection
			allSubprobs[s] = oneProblem()
			allSubprobs[s] = sp_formulation(power, param, stoc)
		 	allSubprobs[s] = attach_scenario(allSubprobs[s], stoc, [s], exargs[:MODEL], 0.0, exargs, subprobType="tight")
		end
		for s in selection
			oneCheckFeaTime = @elapsed isolate_check_one_design_feasibility(power, param, stoc, exargs, stoc.sbdColumns[s], stoc.S, builtModel = allSubprobs)
			isoDesignPool[s] = stoc.sbdColumns[s]
			info("[ISO] Scen $s -> cost = [$(round.(isoDesignPool[s].cost,2))][LB=$(round.(isoDesignPool[s].lb,2))][TIME=$(round.(isoDesignPool[s].time,2))][Cover $(round.(isoDesignPool[s].coverage,2))][CHECKFEATime $(round.(oneCheckFeaTime,2))]")
			checkFeaTime += oneCheckFeaTime
		end
	end

	isoTime = solveTime + checkFeaTime
	info("[TICTOC] main process over in $(isoTime)s")

	for s in selection
		if isoDesignPool[s].cost <= earlyExitCost
			earlyExitFeaCnt = sum(isoDesignPool[s].feamap)
			earlyExitCost = isoDesignPool[s].cost
			earlyExitIndex = s
		end
		push!(stoc.scenarios[s].pool, isoDesignPool[s])
		isoDesignCover[s] = isoDesignPool[s].coverage
		isoDesignCost[s] = isoDesignPool[s].cost
		isoDesignLB[s] = isoDesignPool[s].lb
	end

	# Early exit rule :: if minimum cost isolated scenario fits the risk constraint, then find feasible lower bound.
	if earlyExitFeaCnt >= S * (1-exargs[:eps])
		info("[ISO] Early stopping at isolating stage. Resulting scenario = $earlyExitIndex")
		info("[ISO] Early stopping objective = $earlyExitCost")
		earlyExit = true
	end

	# ============================== Security stage =============================== #
	# Notice this security stage might violate certain first stage solutions
	info("[ISO] Pumping a safety column into decision space")
	isoUnionDesign = union_design(isoDesignPool, param=param)
	isoUnionDesign.feamap = ones(Int, S)
	isoUnionDesign.source = [1:S;]
	push!(stoc.sbdColumns, isoUnionDesign)
	# ============================================================================== #

	return isoTime, stoc, isoUnionDesign.cost, isoDesignCost, earlyExit
end

"""
	Modular Function: Solve one scenario | Fits parallel structure
"""
function isolate_solve_one_scenario(power::Dict,
									param::Dict,
									stoc::stocType,
									s::Int,
									exargs::Dict,
									formulation::Function,
									subprobType::AbstractString;
									kwargs...)

	oneIsoProb = formulation(power,param,stoc, [s], exargs, subprobType=subprobType)
	warmstart_heuristic(oneIsoProb, power, param, stoc, exargs, selection=[s])
	config_solver(oneIsoProb.model, license=config.ENVS, timelimit=config.TIMELIMITIII, mipgap=config.OPTGAP, showlog=config.SHOWLOG, focus="optimality", presolve=1, threads=config.WORKERTHREADS)
	status = solve(oneIsoProb.model)

	if status == :Infeasible
		info("[ISO] Creating infea column on scenario $s")
		oneIsoDesign = infea_design(s, param)
	else
		oneIsoDesign = get_design(oneIsoProb.model)
		oneIsoDesign.source = [s]
		oneIsoDesign.time = getsolvetime(oneIsoProb.model)
		oneIsoDesign.lb = oneIsoProb.model.objBound
		oneIsoDesign.active = true
	end

	return oneIsoDesign
end

function isolate_check_one_design_feasibility(power::Dict,
												param::Dict,
												stoc::stocType,
												exargs::Dict,
												oneIsoDesign::designType,
												S::Int;
												kwargs...)

	options = Dict(kwargs)
	if haskey(options, :builtModel)
		non, oneIsoDesign.coverage, non, feaPool = check_feasible(power, param, stoc, oneIsoDesign, exargs, [], precheck=true, builtModel=options[:builtModel])
	else
		non, oneIsoDesign.coverage, non, feaPool = check_feasible(power, param, stoc, oneIsoDesign, exargs, [], precheck=true)
	end

	oneIsoDesign.feamap = zeros(Int, S)
	for subs in feaPool
		oneIsoDesign.feamap[subs] = 1
	end

	return oneIsoDesign
end
