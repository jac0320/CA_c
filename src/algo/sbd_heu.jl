function improver_heu(power::Dict, param::Dict, stoc::stocType, exargs::Dict,
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
		println("[CG-MIP] Fetching parameter setup from improver_heuristic...")
		if exargs[:CGMAX] > 0
			return Dict(:noimprovestop=>exargs[:CGMAX])
		else
			return Dict(:noimprovestop=>3)
		end
	end

	println("[CG-IMP]Running CG heuristic :: improver_heuristic")

	fixedLB = extras[:fixedLB]
	S = exargs[:S]

	haskey(extras, :tictoc) ? tictoc = extras[:tictoc] : tictoc = 0.0  	#Restart the clock for nothing
	haskey(extras, :incumbent) ? incumbent = extras[:incumbent] : incumbent = inf

	(isa(kickSet, Array)) && (blendSet = collect(blendSet))
	(isa(blendSet, Array)) && (blendSet = collect(blendSet))

	if config.PARALLEL

		# ====================== Spinning Scenario Phase ======================= #
		println("[SPINNING]pickScenarioPool is $pickScenarioPool")
		spinningPairPool = lineup_generate_pool(stoc, pickScenarioPool, pickScenarioPool, solved)

		for pair in spinningPairPool
			solved[pair] = 1
		end

		# Prepare a solution space for carrying design solutions
		spinningPairDesigns = Array{designType}(length(spinningPairPool))

		# Solve each pair in parallel, utilize as many as workers available
		solvePairTime=@elapsed spinningPairDesigns = pmap((a1,a2,a3,a4,a5,a6)->lineup_solve_pair(a1,a2,a3,a4,a5,a6),
														[power for pair in spinningPairPool],
														[param for pair in spinningPairPool],
														[stoc for pair in spinningPairPool],
														[exargs for pair in spinningPairPool],
														[pair for pair in spinningPairPool],
														[subprob_formulation for pair in spinningPairPool])
		println("[SPINNING] Solving paried pool took $(solvePairTime)")
		tictoc += solvePairTime
		println("[TICTOC] Finished spinning solving phase TIME = $tictoc")

		# Regroup the calculated designs for feasibility check
		spinningGroupPairDesigns = lineup_regroup_design(spinningPairDesigns)
		println("[SPINNING] Paired design pool was parititioned into $(length(spinningGroupPairDesigns)) groups")
		feaCheckTime = @elapsed checkedDesigns = pmap((a1,a2,a3,a4,a5)->lineup_groupdesign_fea(a1,a2,a3,a4,a5),
														[power for group in spinningGroupPairDesigns],
														[param for group in spinningGroupPairDesigns],
														[stoc for group in spinningGroupPairDesigns],
														[exargs for group in spinningGroupPairDesigns],
														[group for group in spinningGroupPairDesigns])
		println("[SPINNING] Checking generated design's feasibility took $(feaCheckTime)")
		tictoc += feaCheckTime
		println("[TICTOC] After finishing spinning post process TIME = $tictoc")

		for g in checkedDesigns
			for d in g
				println("([Spinning-S] $(d.source) cost = <$(round.(d.cost,2))($(round.(d.coverage*100,2))\%)> [LB >> $(round.(d.lb,2))];\n")
				stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
			end
		end

		# ====================== Blending Scenarios Phase ====================== #
		println("[BLENDING] blendSet is $blendSet")
		# Generate a set of pairs that are required to be solved
		blendPairPool = lineup_generate_pool(stoc, blendSet, pickScenarioPool, solved)
		for pair in blendPairPool
			solved[pair] = 1
		end

		# Prepare a solution space for carrying design solutions
		blendPairDesigns = Array{designType}(length(blendPairPool))

		# Solve each pair in parallel, utilize as many as workers available
		solvePairTime=@elapsed blendPairDesigns = pmap((a1,a2,a3,a4,a5,a6)->lineup_solve_pair(a1,a2,a3,a4,a5,a6),
														[power for pair in blendPairPool],
														[param for pair in blendPairPool],
														[stoc for pair in blendPairPool],
														[exargs for pair in blendPairPool],
														[pair for pair in blendPairPool],
														[subprob_formulation for pair in blendPairPool])

		println("[BLENDING] Solving paried pool took $(solvePairTime)")
		tictoc += solvePairTime
		println("[TICTOC] After finished blending solving phase TIME = $tictoc")

		# Regroup the calculated designs for feasibility check
		blendGroupPairDesigns = lineup_regroup_design(blendPairDesigns)
		println("[BLENDING] Paired design pool was parititioned into $(length(blendGroupPairDesigns)) groups")

		# Conduct feasibility check in parallel, warm start scheme can be applied here
		feaCheckTime=@elapsed blendGroupPairDesigns = pmap((a1,a2,a3,a4,a5)->lineup_groupdesign_fea(a1,a2,a3,a4,a5),
														[power for group in blendGroupPairDesigns],
														[param for group in blendGroupPairDesigns],
														[stoc for group in blendGroupPairDesigns],
														[exargs for group in blendGroupPairDesigns],
														[group for group in blendGroupPairDesigns])
		println("[BLENDING] Checking generated design's feasibility took $(feaCheckTime)")
		tictoc += feaCheckTime
		println("[TICTOC] After finishing blending post process TIME = $tictoc")

		# Collect all accumulated columns into the public pool
		for g in blendGroupPairDesigns
			for d in g
				println("($(myid()))[Blending-S] $(d.source) cost = <$(round.(d.cost,2))($(round.(d.coverage*100,2))\%)> [LB >> $(round.(d.lb,2))];\n")
				stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
			end
		end

	else
		# =============================== Sequential Mode =============================== #

		# ====================== Spinning Scenario Phase ======================= #
		println("[SPINNING]pickScenarioPool is $pickScenarioPool")
		# Generate a set of pairs that are required to be solved
		spinningPairPool = lineup_generate_pool(stoc, pickScenarioPool, pickScenarioPool, solved)

		println("SPINNING TODO Counts: $(length(spinningPairPool))")
		# Prepare a solution space for carrying design solutions
		spinningPairDesigns = Array{designType}(length(spinningPairPool))

		# Solve each pair sequentially and store the design in an array
		idx = 0
		solveTime = 0
		for pair in spinningPairPool
			idx += 1
			idxSolveTime = @elapsed spinningPairDesigns[idx] = lineup_solve_pair(power, param, stoc, exargs, pair, subprob_formulation)
			solveTime += idxSolveTime
			solved[pair] = 1
		end
		println("[SPINNING] Solving paried pool took $(solveTime)")

		# Check each generated design sequentially
		feaCheckTime = 0.0
		if config.WARMSTART
			allSubprobs = Array{oneProblem}(stoc.S)
			# Prepare a vector of problems for repetitvely checking
			for ss = 1:stoc.S
				allSubprobs[ss] = oneProblem()
				allSubprobs[ss] = sbd_base_formulation(power, param, stoc)
				allSubprobs[ss] = attach_scenario(allSubprobs[ss], stoc, [ss], exargs[:MODEL], 0.0, exargs, subprobType="tight")
			end
			for i in 1:length(spinningPairDesigns)
				singleCheckTime = @elapsed spinningPairDesigns[i] = isolate_check_one_design_feasibility(power, param, stoc, exargs, spinningPairDesigns[i], stoc.S, builtModel = allSubprobs)
				d = spinningPairDesigns[i]
				println("($(myid()))[Spinning-S] $(d.source) cost = <$(round.(d.cost,2))(Cover $(round.(d.coverage,2)))(Time $(round.(d.time,2))s\%)> [LB >> $(round.(d.lb,2))];")
				feaCheckTime += singleCheckTime
			end
		else
			error("Should always have problem built.")
			for i in 1:length(spinningPairDesigns)
				singleCheckTime = @elapsed spinningPairDesigns[i] = isolate_check_one_design_feasibility(power,param, stoc,exargs,spinningPairDesigns[i], stoc.S)
				feaCheckTime += singleCheckTime
			end
		end
		println("[SPINNING] Checking generated design's feasibility took $(feaCheckTime)")

		# Collect all accumulated columns into the public pool
		for d in spinningPairDesigns
			# println("($(myid()))[Spinning-S] $(d.source) cost = <$(round.(d.cost,2))($(round.(d.time,2))s\%)> [LB >> $(round.(d.lb,2))];")
			stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
		end

		# ====================== Blending Scenario Phase ======================= #
		println("[BLENDING]blendSet is $blendSet")
		# Generate a set of pairs that are required to be solved
		blendPairPool = lineup_generate_pool(stoc, blendSet, pickScenarioPool, solved)
		println("[BLENDING] Blending phase TODO counts $(length(blendPairPool))")
		# Prepare a solution space for carrying design solutions
		blendPairDesigns = Array{designType}(length(blendPairPool))

		# Solve each pair sequentially and store the design in an array
		idx = 0
		solveTime = 0
		for pair in blendPairPool
			idx += 1
			idxSolveTime = @elapsed blendPairDesigns[idx] = lineup_solve_pair(power, param, stoc, exargs, pair, subprob_formulation)
			solveTime += idxSolveTime
			solved[pair] = 1
		end
		println("[BLENDING] Solving paried pool took $(solveTime)")

		# Check each generated design sequentially
		feaCheckTime = 0.0
		if config.WARMSTART
			# allSubprobs already created
			for i in 1:length(blendPairDesigns)
				singleCheckTime = @elapsed blendPairDesigns[i] = isolate_check_one_design_feasibility(power, param, stoc, exargs, blendPairDesigns[i], stoc.S, builtModel = allSubprobs)
				feaCheckTime += singleCheckTime
			end
		# else
		# 	for i in 1:length(blendPairDesigns)
		# 		singleCheckTime = @elapsed blendPairDesigns[i] = isolate_check_one_design_feasibility(power,param, stoc,exargs,blendPairDesigns[i], stoc.S)
		# 		feaCheckTime += singleCheckTime
		# 	end
		end
		println("[BLENDING] Checking generated design's feasibility took $(feaCheckTime)")

		# Collect all accumulated columns into the public pool
		for d in blendPairDesigns
			println("($(myid()))[BLENDING-S] $(d.source) cost = <$(round.(d.cost,2))($(round.(d.coverage*100,2))\%)> [LB >> $(round.(d.lb,2))];")
			stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
		end

	end

	return stoc, solved, tictoc
end

function partial_heu(power::Dict, param::Dict, stoc::stocType, exargs::Dict,
									masterIncumbent::Float64,
									kickSet,
									blendSet,
									pickScenarioPool::Array,
									neglectedScenarioPool::Array,
									master_formulation::Function,
									subprob_formulation::Function,
									solved::Dict, allSubprobs=nothing;
									kwargs...)

	extras = Dict(kwargs)

	if haskey(extras, :findparameter)
		println("Fetching parameter setup from improver_heuristic...")
		if exargs[:CGMAX] > 0
			return Dict(:noimprovestop=>exargs[:CGMAX])
		else
			return Dict(:noimprovestop=>6)
		end
	end

	println("Running CG heuristic :: lineup_partial_heuristic")
	fixedLB = extras[:fixedLB]
	S = exargs[:S]

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
		println("[PARTIAL]pickScenarioPool is $pickScenarioPool")

		# Collect a set of iso costs
		isoCosts = []
		isoIdx = [] # The mismatch resulted from master problem selection
		println("[PARTIAL] CG Heuristic response on iteration $iterIdxer")
		restIdx = []
		restCosts = []
		for s in 1:exargs[:S]
			if s in pickScenarioPool
				push!(isoCosts, stoc.sbdColumns[s].cost)
				push!(isoIdx, s)
			else
				push!(restCosts, stoc.sbdColumns[s].cost)
				push!(restIdx, s)
			end
		end
		println("[PARTIAL] CG Heuristic collect isoCosts $isoCosts")
		println("[PARTIAL] CG Heuristic coressponding idx $isoIdx")

		if iterIdxer > length(isoCosts)
			rankedScenario = restIdx[findfirst(restCosts, select(restCosts, iterIdxer-length(isoCosts)))]
			println("[PARTIAL]Picking scenario $rankedScenario for this iteration, outside the alpha scope")
		else
			rankedScenario = isoIdx[findfirst(isoCosts, select(isoCosts, iterIdxer))]
			println("[PARTIAL]Picking scenario $rankedScenario for this iteration")
		end

		# Generate a set of pairs that are required to be solved
		partialPairPool = partial_generate_pool(stoc, rankedScenario, [pickScenarioPool;blendSet], solved)
		println("[PARTIAL] TODO List: $partialPairPool")

		# Extra synced record for parallel scheme
		for pair in partialPairPool
			solved[pair] = 1
		end

		# Prepare a solution space for carrying design solutions
		partialPairDesigns = Array{designType}(length(partialPairPool))

		# Solve each pair in parallel, utilize as many as workers available
		solvePairTime=@elapsed partialPairDesigns = pmap((a1,a2,a3,a4,a5,a6)->lineup_solve_pair(a1,a2,a3,a4,a5,a6),
														[power for pair in partialPairPool],
														[param for pair in partialPairPool],
														[stoc for pair in partialPairPool],
														[exargs for pair in partialPairPool],
														[pair for pair in partialPairPool],
														[subprob_formulation for pair in partialPairPool])

		# Finish Solving
		println("[PARTIAL] Solving paried pool took $(solvePairTime)")
		tictoc += solvePairTime
		println("[PARTIAL] Finished partial solving phase TIME = $tictoc")

		# Regroup the calculated designs for feasibility check
		partialGroupPairDesigns = lineup_regroup_design(partialPairDesigns)
		println("[SHORT] Paired design pool was parititioned into $(length(partialGroupPairDesigns)) groups")
		feaCheckTime = @elapsed checkedDesigns = pmap((a1,a2,a3,a4,a5)->lineup_groupdesign_fea(a1,a2,a3,a4,a5),
														[power for group in partialGroupPairDesigns],
														[param for group in partialGroupPairDesigns],
														[stoc for group in partialGroupPairDesigns],
														[exargs for group in partialGroupPairDesigns],
														[group for group in partialGroupPairDesigns])
		println("[PARTIAL] Checking generated design's feasibility took $(feaCheckTime)")
		tictoc += feaCheckTime
		println("[TICTOC] After finishing partial post process TIME = $tictoc")

		# Collect all accumulated columns into the public pool
		for g in checkedDesigns
			for d in g
				println("($(myid()))[PARTIAL-S] $(d.source) cost = <$(round.(d.cost,2))($(round.(d.coverage*100,2))\%)> [LB >> $(round.(d.lb,2))];\n")
				stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
			end
		end

	else
		# =============================== Sequential Mode =============================== #

		# ====================== Spinning Scenario Phase ======================= #
		println("[PARTIAL]Master picked $pickScenarioPool")

		# Collect a set of iso costs
		isoCosts = []
		isoIdx = [] # The mismatch resulted from master problem selection
		println("[PARTIAL] CG Heuristic response on iteration $iterIdxer")
		restIdx = []
		restCosts = []
		for s in 1:exargs[:S]
			if s in pickScenarioPool
				push!(isoCosts, stoc.sbdColumns[s].cost)
				push!(isoIdx, s)
			else
				push!(restCosts, stoc.sbdColumns[s].cost)
				push!(restIdx, s)
			end
		end
		println("[PARTIAL] CG Heuristic collect isoCosts $isoCosts")
		println("[PARTIAL] CG Heuristic coressponding idx $isoIdx")

		if iterIdxer > length(isoCosts)
			rankedScenario = restIdx[findfirst(restCosts, select(restCosts, iterIdxer-length(isoCosts)))]
			println("[PARTIAL]Picking scenario $rankedScenario for this iteration, outside the alpha scope")
		else
			rankedScenario = isoIdx[findfirst(isoCosts, select(isoCosts, iterIdxer))]
			println("[PARTIAL]Picking scenario $rankedScenario for this iteration")
		end

		# Generate a set of pairs that are required to be solved
		# partialPairPool = partial_generate_pool(stoc, pickScenarioPool[iterIdxer%(length(pickScenarioPool))+1], pickScenarioPool, solved)
		partialPairPool = partial_generate_pool(stoc, rankedScenario, [pickScenarioPool;blendSet], solved)

		println("[PARTIAL] TODO List: $partialPairPool")
		partialPairDesigns = Array{designType}(length(partialPairPool))

		# Solve each pair sequentially and store the design in an array
		idx = 0
		solveTime = 0
		for pair in partialPairPool
			idx += 1
			idxSolveTime = @elapsed partialPairDesigns[idx] = lineup_solve_pair(power, param, stoc, exargs, pair, subprob_formulation)
			solveTime += idxSolveTime
			solved[pair] = 1
		end
		println("[PARTIAL] Solving paried pool took $(solveTime)")

		# Check each generated design sequentially
		feaCheckTime = 0.0
		if config.WARMSTART
			# allSubprobs = Array{oneProblem}(stoc.S)
			# # Prepare a vector of problems for repetitvely checking
			# for ss = 1:stoc.S
			# 	allSubprobs[ss] = oneProblem()
			# 	allSubprobs[ss] = sbd_base_formulation(power, param, stoc)
			# 	allSubprobs[ss] = attach_scenario(allSubprobs[ss], stoc, [ss], exargs[:MODEL], 0.0, exargs, subprobType="tight")
			# end
			for i in 1:length(partialPairDesigns)
				singleCheckTime = @elapsed partialPairDesigns[i] = isolate_check_one_design_feasibility(power, param, stoc, exargs, partialPairDesigns[i], stoc.S, builtModel = allSubprobs)
				d = partialPairDesigns[i]
				println("($(myid()))[Spinning-S] $(d.source) cost = <$(round.(d.cost,2))(Cover $(round.(d.coverage,2)))(Time $(round.(d.time,2))s\%)> [LB >> $(round.(d.lb,2))];")
				feaCheckTime += singleCheckTime
			end
		else
			for i in 1:length(partialPairDesigns)
				singleCheckTime = @elapsed partialPairDesigns[i] = isolate_check_one_design_feasibility(power,param, stoc,exargs,partialPairDesigns[i], stoc.S)
				feaCheckTime += singleCheckTime
			end
		end
		println("[PARTIAL] Checking generated design's feasibility took $(feaCheckTime)")

		# Collect all accumulated columns into the public pool
		for d in partialPairDesigns
			stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, d)
		end

	end

	return stoc, solved, tictoc
end
