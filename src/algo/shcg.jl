function shcg(param::Dict, stoc::stocType, driver::Dict, master_formulation::Function, subprob_formulation::Function)

	# Algorithm Initialization
	S = param[:S]
	eps = param[:eps]
	st = time()
	iter = 0
	noImprove = 0
	incumbub = Inf
	incumblb = Inf
	incumbgap = Inf
	incumbcol = nothing
	solved = Dict()
	cgparameters = eval(parse(driver[:CGHEURISTIC]))(driver, findparameter=true)
	isolate_stage(param, stoc, driver)
	isoCost = [stoc.sbdColumns[s].cost for s in 1:S]
	incumblb = select(isoCost, convert(Int, floor((1-eps)*S)))
	ubpool = [s for s in 1:S if isoCost[s] < incumblb]
	heuristic_col = shcg_solve_sp(param, stoc, driver, ubpool, isoCost)
	shcg_store_col(stoc.sbdColumns, heuristic_col)
	incumbcol = search_incumb_col(stoc.sbdColumns, param)
	incumbub = incumbcol.cost
	builtmodel = [build_sp(param, stoc, driver, selection=[s], sbtype="tight") for s in 1:S]

	while true

		iter += 1 # Update iterator

		info("[SBD] >>>>>>> Iteration $iter :: $incumbent")
		columns_cleaner(stoc, incumbent)

		P = length([i for i in 1:length(stoc.sbdColumns) if stoc.sbdColumns[i].active])
		info("[SBD] Active Count = $(P)")

		mastersol = shcg_solve_master(param, stoc, driver, prevsol=incumbcol, incumb=mastersol, builtmodel=builtmodel)
		candidatecol = search_incumb_col(stoc.sbdColumns, param)

		if candidatecol.cost < incumbub
			incumbcol = candidatecol
			incumbub = candidatecol.cost
			incumbgap = (incumbub - incumblb)/incumbub
			info("[SBD] Improvement detected. Current incumbent $incumbent[$incumbentGap][Gap $(round.((incumbent-incumblb)/incumbent*100,2))%]")
			noImprove = 0
		else
			noImprove += 1
			info("[SBD] NO Improvement.Current incumbent $incumbent[$incumbentGap][Gap $(round.((incumbent-incumblb)/incumbent*100,2))%]")
		end

		# Termination Conditions
		signal = :Continue
		if noImprove > cgparameters[:noimprovestop]
			info("[SHCG] TERMINATION : No improvement strike out")
			signal = :UserLimits
		elseif iter > driver[:S]*(1-driver[:EPS])
			info("[SHCG] TERMINATION : Maximum Iteartion")
			signal = :UserLimits
		elseif time()-st > dirver[:TIMELIMIT]
			info("[SHCG] TERMINATION : Timeout")
			signal = :UserLimits
		elseif length([i for i in 1:S if isoCost[i] > incumbub]) >= ceil(eps*S)
			info("[SHCG] TERMINATION : Enough scenario eliminated")
			signal = :Optimal
			postoptpool = [i for i in 1:S if isoCost[i] <= incumbub]
		elseif incumbgap <= driver[:OPTGAP]
			info("[SHCG] TERMINATION : Gap closed")
			signal = :Optimal
			postoptpool = [i for i in 1:S if incumbcol.feamap[i] == 1]
		end

		# Terminate the algorithm
		if signal != :Continue
			info("[POST] Starting POST optimization process")

			postoptcol = shcg_solve_sp(param, stoc, driver, postoptpool, isoCost)
			shcg_store_col(stoc.sbdColumns, postoptcol)

			outputcol = search_incumb_col(stoc.sbdColumns, param)
			outputlb = incumblb
			outputub = outputcol.cost
			outputgap = (outputub - outputlb)/outputub
			print_design(outputcol, driver)
			write_output_files(sbdDesign, driver)
			info("[EXIT] COL-IDX $(outputcol.k) || UB $(outputub) || LB $(outputlb) || GAP $(outputgap)")
			info("[EXIT] Finish time $(time()-st) s")
			info("[EXIT] Successfully completed SBD heuristic.")

			return outputcol
		end

		eval(parse(driver[:CGHEURISTIC]))(param, stoc, driver,
										  solved=solved,
										  incumb=incumbub,
										  iter=iter,
										  builtmodel=builtmodel)
	end
end

function shcg_enu(driver::Dict; findparameter=false)

	if findparameter
		if driver[:CGMAX] > 0
			return Dict(:noimprovestop=>driver[:CGMAX])
		else
			return Dict(:noimprovestop=>3)
		end
	end

	return Dict()
end

function shcg_enu(param::Dict, stoc::stocType, driver::Dict; kwargs...)

	options = Dict(kwargs)

	solved = options[:solved]
	incumb = options[:incumb]
	iter = options[:iter]

	jobs = shcg_generate_task(param, solved, dim=iter)
	info("[ENU] Job Count = $(length(jobs))")

	if driver[:PARALLEL]
		tSolve = @elapsed cols = pmap((a1,a2,a3,a4)->shcg_solve_sp(a1,a2,a3,a4),
										[param for j in jobs],
										[stoc for j in jobs],
										[driver for j in jobs],
										[j for j in jobs])
		info("[ENU] Solve time $(tsolve)")

		groups = shcg_regroup_cols(cols, driver)
		tcheck = @elapsed groups = pmap((a1,a2,a3,a4)->shcg_check_fea(a1,a2,a3,a4),
										[param for g in groups],
										[stoc for g in groups],
										[driver for g in groups],
										[g for g in groups])
		info("[ENU] Feasibility check time $(tcheck)")

		tstore = @elapsed [shcg_store_col(stoc.sbdColumns, g, incumb=incumb) for g in groups]
		info("[ENU] Storaging time $(tstore)")s
	else
		tsolve = @elapsed cols = [shcg_solve_sp(param, stoc, driver, j) for j in jobs]
		info("[ENU] Solve time $(tsolve)")

		feamodels = [build_sp(param, stoc, driver, selection=s, sbtype="tight") for s in subsets([1:S;], 1)]
		info("[ENU] Built $(length(feamodels)) feaibility check models.")

		tcheck = @elapsed cols = [iso_check_block(c, param, stoc, driver, builtmodel=feamodels) for c in cols]
		info("[ENU] Feasibility check time $(tcheck)")

		tstore = @elapsed [shcg_store_col(stoc.sbdColumns, c, incumb=incumb) for c in cols]
		info("[ENU] Storaging time $(tstore)")
	end

	return
end

function partial_heu(power::Dict, param::Dict, stoc::stocType, driver::Dict,
									masterIncumbent::Float64,
									kickSet,
									blendSet,
									pickScenarioPool::Array,
									neglectedScenarioPool::Array,
									master_formulation::Function,
									subprob_formulation::Function,
									solved::Dict;
									kwargs...)

	extras = Dict(kwargs)

	if haskey(extras, :findparameter)
		info("Fetching parameter setup from improver_heuristic...")
		if driver[:CGMAX] > 0
			return Dict(:noimprovestop=>driver[:CGMAX])
		else
			return Dict(:noimprovestop=>6)
		end
	end

	info("Running CG heuristic :: lineup_partial_heuristic")
	incumblb = extras[:incumblb]
	S = driver[:S]

	if haskey(extras, :tictoc)
		tictoc = extras[:tictoc]
	else
		tictoc = 0.0  	#Restart the clock for nothing
 	end

	if haskey(extras, :incumbent)
		incumbent = extras[:incumbent]
	else
		incumbent = inf
	end

	if haskey(extras, :iterIdxer)
		iterIdxer = extras[:iterIdxer]
	end

	# Change of data structure if not array
	if !isa(kickSet, Array)
		blendSet = collect(blendSet)
	end

	if !isa(blendSet, Array)
		blendSet = collect(blendSet)
	end

	# Parallel Implementation
	if config.PARALLEL

		# ====================== Spinning Scenario Phase ======================= #
		info("[PARTIAL]pickScenarioPool is $pickScenarioPool")

		# Collect a set of iso costs
		isoCosts = []
		isoIdx = [] # The mismatch resulted from master problem selection
		info("[PARTIAL] CG Heuristic response on iteration $iterIdxer")
		restIdx = []
		restCosts = []
		for s in 1:driver[:S]
			if s in pickScenarioPool
				push!(isoCosts, stoc.sbdColumns[s].cost)
				push!(isoIdx, s)
			else
				push!(restCosts, stoc.sbdColumns[s].cost)
				push!(restIdx, s)
			end
		end
		info("[PARTIAL] CG Heuristic collect isoCosts $isoCosts")
		info("[PARTIAL] CG Heuristic coressponding idx $isoIdx")

		if iterIdxer > length(isoCosts)
			rankedScenario = restIdx[findfirst(restCosts, select(restCosts, iterIdxer-length(isoCosts)))]
			info("[PARTIAL]Picking scenario $rankedScenario for this iteration, outside the alpha scope")
		else
			rankedScenario = isoIdx[findfirst(isoCosts, select(isoCosts, iterIdxer))]
			info("[PARTIAL]Picking scenario $rankedScenario for this iteration")
		end

		# Generate a set of pairs that are required to be solved
		partialPairPool = partial_generate_pool(stoc, rankedScenario, [pickScenarioPool;blendSet], solved)
		info("[PARTIAL] TODO List: $partialPairPool")

		# Extra synced record for parallel scheme
		for pair in partialPairPool
			solved[pair] = 1
		end

		# Prepare a solution space for carrying design solutions
		partialPairDesigns = Array{designType}(length(partialPairPool))

		# Solve each pair in parallel, utilize as many as workers available
		solvePairTime=@elapsed partialPairDesigns = pmap((a1,a2,a3,a4,a5,a6)->shcg_solve_sp(a1,a2,a3,a4,a5,a6),
														[power for pair in partialPairPool],
														[param for pair in partialPairPool],
														[stoc for pair in partialPairPool],
														[driver for pair in partialPairPool],
														[pair for pair in partialPairPool],
														[subprob_formulation for pair in partialPairPool])

		# Finish Solving
		info("[PARTIAL] Solving paried pool took $(solvePairTime)")
		tictoc += solvePairTime
		info("[PARTIAL] Finished partial solving phase TIME = $tictoc")

		# Regroup the calculated designs for feasibility check
		partialGroupPairDesigns = lineup_regroup_design(partialPairDesigns)
		info("[SHORT] Paired design pool was parititioned into $(length(partialGroupPairDesigns)) groups")
		feaCheckTime = @elapsed checkedDesigns = pmap((a1,a2,a3,a4,a5)->shcg_check_fea(a1,a2,a3,a4,a5),
														[power for group in partialGroupPairDesigns],
														[param for group in partialGroupPairDesigns],
														[stoc for group in partialGroupPairDesigns],
														[driver for group in partialGroupPairDesigns],
														[group for group in partialGroupPairDesigns])
		info("[PARTIAL] Checking generated design's feasibility took $(feaCheckTime)")
		tictoc += feaCheckTime
		info("[TICTOC] After finishing partial post process TIME = $tictoc")

		# Collect all accumulated columns into the public pool
		for g in checkedDesigns
			for d in g
				info("($(myid()))[PARTIAL-S] $(d.source) cost = <$(round.(d.cost,2))($(round.(d.coverage*100,2))\%)> [LB >> $(round.(d.lb,2))];\n")
				stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
			end
		end

	else
		# =============================== Sequential Mode =============================== #

		# ====================== Spinning Scenario Phase ======================= #
		info("[PARTIAL]Master picked $pickScenarioPool")

		# Collect a set of iso costs
		isoCosts = []
		isoIdx = [] # The mismatch resulted from master problem selection
		info("[PARTIAL] CG Heuristic response on iteration $iterIdxer")
		restIdx = []
		restCosts = []
		for s in 1:driver[:S]
			if s in pickScenarioPool
				push!(isoCosts, stoc.sbdColumns[s].cost)
				push!(isoIdx, s)
			else
				push!(restCosts, stoc.sbdColumns[s].cost)
				push!(restIdx, s)
			end
		end
		info("[PARTIAL] CG Heuristic collect isoCosts $isoCosts")
		info("[PARTIAL] CG Heuristic coressponding idx $isoIdx")

		if iterIdxer > length(isoCosts)
			rankedScenario = restIdx[findfirst(restCosts, select(restCosts, iterIdxer-length(isoCosts)))]
			info("[PARTIAL]Picking scenario $rankedScenario for this iteration, outside the alpha scope")
		else
			rankedScenario = isoIdx[findfirst(isoCosts, select(isoCosts, iterIdxer))]
			info("[PARTIAL]Picking scenario $rankedScenario for this iteration")
		end

		# Generate a set of pairs that are required to be solved
		# partialPairPool = partial_generate_pool(stoc, pickScenarioPool[iterIdxer%(length(pickScenarioPool))+1], pickScenarioPool, solved)
		partialPairPool = partial_generate_pool(stoc, rankedScenario, [pickScenarioPool;blendSet], solved)

		info("[PARTIAL] TODO List: $partialPairPool")
		partialPairDesigns = Array{designType}(length(partialPairPool))

		# Solve each pair sequentially and store the design in an array
		idx = 0
		solveTime = 0
		for pair in partialPairPool
			idx += 1
			idxSolveTime = @elapsed partialPairDesigns[idx] = shcg_solve_sp(power, param, stoc, driver, pair, subprob_formulation)
			solveTime += idxSolveTime
			solved[pair] = 1
		end
		info("[PARTIAL] Solving paried pool took $(solveTime)")

		# Check each generated design sequentially
		feaCheckTime = 0.0
		if config.WARMSTART
			allSubprobs = Array{oneProblem}(stoc.S)
			# Prepare a vector of problems for repetitvely checking
			for ss = 1:stoc.S
				allSubprobs[ss] = oneProblem()
				allSubprobs[ss] = sp_formulation(power, param, stoc)
				allSubprobs[ss] = attach_scenario(allSubprobs[ss], stoc, [ss], driver[:MODEL], 0.0, driver, subprobType="tight")
			end
			for i in 1:length(partialPairDesigns)
				singleCheckTime = @elapsed partialPairDesigns[i] = iso_check_block(power, param, stoc, driver, partialPairDesigns[i], stoc.S, builtmodel = allSubprobs)
				d = partialPairDesigns[i]
				info("($(myid()))[Spinning-S] $(d.source) cost = <$(round.(d.cost,2))(Cover $(round.(d.coverage,2)))(Time $(round.(d.time,2))s\%)> [LB >> $(round.(d.lb,2))];")
				feaCheckTime += singleCheckTime
			end
		else
			for i in 1:length(partialPairDesigns)
				singleCheckTime = @elapsed partialPairDesigns[i] = iso_check_block(power,param, stoc,driver,partialPairDesigns[i], stoc.S)
				feaCheckTime += singleCheckTime
			end
		end
		info("[PARTIAL] Checking generated design's feasibility took $(feaCheckTime)")

		# Collect all accumulated columns into the public pool
		for d in partialPairDesigns
			stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
		end

	end

	return stoc, solved, tictoc
end

function shcgnr(param::Dict, stoc::stocType, driver::Dict; selection=[], isoCost=[])

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	S, allS = length(selection), param[:S]

	iter = 0
	ast = time()
	scenpool = Int[]

	if isempty(isoCost)
		st = time()
		isolate_stage(param, stoc, driver, selection=selection, skipsafety=true)
		info("[SBD-NR] ISO stage $(time()-st) || A-T $(time()-ast)")
	end

	builtmodel = [build_sp(param, stoc, driver, selection=[s], sbtype="slackness") for s in 1:allS]

	isoCost = [i.cost for i in stoc.sbdColumns]
	nextidx = findfirst(isoCost, maximum(isoCost))
	incumbCol = stoc.sbdColumns[nextidx]
	push!(scenpool, nextidx)

	info("[SBD-NR] Init -- $(scenpool)")
	infeas = [s for s in selection if incumbCol.feamap[s] == 0]

	while true
		iter += 1
		st = time()

		nextidx, infeas = check_slackness(param, stoc, incumbCol, driver, selection=infeas, builtmodel=builtmodel)
		info("[SBD-NR] CHECK SLACKNESS $(time()-st) || A-T $(time()-ast)")

		cover = (S - length(infeas))/S*100
		info("[SBD-NR] FEACNT $(S-length(infeas)) || COVERAGE $(round.(cover,3))%")

		if isempty(infeas) || iter > driver[:MAXITER]
			info("[SBD-NR: EXIT] Termination condition reached...")
			if iter == 1
				info("[SBD-NR: EXIT] POST OPT $(scenpool)")
				postprob = build_sp(param, stoc, driver, selection=scenpool, sbtype="tight")
				# warmstart_heuristic(postprob, stoc, driver, selection=scenpool)
				config_solver(postprob, driver, timelimit=driver[:TIMELIMITII],focus="optimality")
				status = solve(postprob.model, suppress_warnings=true)
				incumbCol = get_design(postprob, idx=iter)
			end
			write_output_files(incumbCol, driver)
			print_design(incumbCol, param)
			return incumbCol
		end

		push!(scenpool, nextidx)	# For next master problem

		info("[SBD-NR] MASETER POOL $scenpool")
		st = time()
		masterprob = build_sp(param, stoc, driver, selection=scenpool, sbtype="tight")
		# warmstart_heuristic(masterprob, param, stoc, driver, selection=scenpool)
		config_solver(masterprob.model, driver, timelimit=driver[:TIMELIMITII], focus="optimality")
		status = solve(masterprob.model, suppress_warnings=true)

		if status == :Infasible  # Assumption for climate problem
			print_iis_gurobi(masterprob.model, driver)
			error("[SBD-NR] Master infeasible")
		end

		incumbCol = get_design(masterprob)
		info("[SBD-NR] ITER $(iter) || UB $(incumbCol.cost) || LB $(incumbCol.lb)")
		info("[SBD-NR] MASTER $(time()-st) || A-TIME $(time()-ast)")
	end
end
