function sbd_heuristic(power::Dict, param::Dict, stoc::stocType, exargs::Dict,
							master_formulation, subprob_formulation; kwargs...)

	println("Solving problem using SBD heuristic method... ")
	println("Setting up solver environments...")
	options = Dict(kwargs)

	incumbent = Inf
	incumbentRecord = []
	lbPool = []
	yesImprove = 0
	noImprove = 0
	stop = false
	iter = 0
	earlyExit = false
	S = exargs[:S]
	eps = exargs[:eps]
	solved = Dict()
	nonCuts = Dict()
	kickSet = []
	recordingCacheISO = false
	tictoc = 0.0
	cuts = []

	cgparameters = eval(parse(exargs[:CGHEURISTIC]))(power, param, stoc, exargs,
		0.0,Set(),Set(),[],[], master_formulation, subprob_formulation, Dict(), findparameter=true)
	btt = time()
	allSubprobs = Array{oneProblem}(stoc.S)
	# Prepare a vector of problems for repetitvely checking
	for ss = 1:stoc.S
		allSubprobs[ss] = oneProblem()
		allSubprobs[ss] = sbd_base_formulation(power, param, stoc)
		allSubprobs[ss] = attach_scenario(allSubprobs[ss], stoc, [ss], exargs[:MODEL], 0.0, exargs, subprobType="tight")
	end
	@show time()-btt
	# ====================================== PHASE 1 ====================================== #
	println("[SBD] Isolating scenarios for basic stage...")
	isoTime, stoc, isoUnionCost, isoCost, earlyExit = isolate_stage(power, param, stoc, exargs, subprob_formulation, [], allSubprobs)
	if earlyExit == true
		return 0 # Here need to return something else
	end
	println("ISO Stage took $(isoTime)s")
	tictoc += isoTime
	println("[ISO] After ISO stage TIME = $tictoc s")

	# Getting a fixed lower bound calculated
	fixedLB = select(isoCost, Int(floor(S*(1-eps))))
	lbScenario = -1
	for s in 1:S
		if abs(isoCost[s] - fixedLB) <= 0.1
			lbScenario = s
			println("[SBD] LB defining scenario is $lbScenario")
		end
	end

	if stoc.sbdColumns[lbScenario].coverage >= (1-eps)
		println("[SBD] LB defining scenario indicates optimality (LB == UB) with coverage [ $(stoc.sbdColumns[lbScenario].coverage)]")
	end

	# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	# Solve a larger problem with all scenarios that may provide good upper bound
	# Select initial UB-Defining Pool
	ubPool = []
	ubISOCost = []
	for s in 1:S
		if isoCost[s] <= fixedLB
			push!(ubPool, s)
			push!(ubISOCost, isoCost[s])
		end
	end
	ubPoolTime = @elapsed terminate, lbPoolProb = sbd_ubPool_solve(power, param, stoc, exargs,
																	fixedLB, ubPool, ubISOCost, subprob_formulation)
	if terminate
		return lbPoolProb, stoc.sbdColumns[end]
	end
	tictoc += ubPoolTime
	println("[TICTOC] After UB Pool TIME = $tictoc")
	# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	# Calculate the initial solution for UB that fits all scenarios
	lbPoolDesign = stoc.sbdColumns[end]
	if lbPoolDesign.cost < isoUnionCost
		incumbent = lbPoolDesign.cost
		incumbentDesign = lbPoolDesign
		incumbentPool = copy(lbPool)
	else
		incumbent = isoUnionCost
		incumbentDesign = stoc.sbdColumns[S+1]
		incumbentPool = [1:S;]
	end

	# ===================== Algorithm main loop ============================== #
	println("[SBD] Entering main loop... INCUMBENT = $incumbent")
	while stop == false # Until termination condition is met

		iter += 1

		println("[SBD] >>>>>>> Iteration $iter :: $incumbent")
		println("[SBD] Cleaning columns before solving master...")
		columns_cleaner(stoc, incumbent)
		totalColumn = 0
		for i = 1:length(stoc.sbdColumns)
			if stoc.sbdColumns[i].active
				totalColumn += 1
			end
		end
		println("[SBD] Accumulated active columns = $totalColumn")

		if iter > 1
			println("[SBD] Solving inital master problem...")
			masterTime = @elapsed masterDesign, pickScenarioPool, neglectedScenarioPool, masterSolution =
				sbd_solve_master(power, param, stoc, exargs, prevSolution=masterSolution, incumbent=incumbent)
		else
			println("[SBD] Solving master problem... with warm start solutions...")
			masterTime = @elapsed masterDesign, pickScenarioPool, neglectedScenarioPool, masterSolution =
				sbd_solve_master(power, param, stoc, exargs)
		end

		gaptofixedLB = (masterDesign.cost - fixedLB)/masterDesign.cost * 100
		println("[SBD] Current gap $gaptofixedLB")

		if masterDesign.cost < incumbent
			incumbent = masterDesign.cost
			incumbentGap = round.(100*(incumbent - fixedLB)/incumbent, 4)
			incumbentPool = copy(pickScenarioPool)
			incumbentDesign = masterDesign
			println("[SBD] Improvement detected. Current incumbent $incumbent[$incumbentGap][Gap $(round.((incumbent-fixedLB)/incumbent*100,2))%]")
			yesImprove += 1
			noImprove = 0
		else
			incumbentGap = round.(100*(incumbent - fixedLB)/incumbent, 4)
			if iter > 2
				noImprove += 1
				println("[SBD] NO Improvement.Current incumbent $incumbent[$incumbentGap][Gap $(round.((incumbent-fixedLB)/incumbent*100,2))%]")
			end
			incumbentPool = copy(incumbentPool)
		end

		write_output_files(power, param, stoc, incumbentDesign, exargs)

		tictoc += masterTime
		println("[TICTOC] After master in iteration $iter, TIME = $tictoc")
		# Termination I :: consecutively no improvement
		stop = false
		if noImprove > cgparameters[:noimprovestop] || iter > exargs[:S]*(1-exargs[:EPS]) || tictoc > 21600
			println("[I] NO INCREMENT/MAX IITERATION STRIKEOUT. Best UB OBTAINED. TERMINATING WITH INUMCBENT = ", incumbent)
			incumbentPool = pickScenarioPool
			stop = true
		end

		# Termination II :: Master objective cuts on the isolated staris
		println("[MAIN] Screening scenarios for inferences... ")
		kickSet = Set()
		blendSet = Set()
		for s in 1:S
			if masterDesign.cost < isoCost[s]
				push!(kickSet, s)
			elseif (s in neglectedScenarioPool)
				push!(blendSet,s)
			end
		end
		kickSet = [collect(kickSet); columns_manager(stoc, kickSet, incumbent, exargs)]
		for s in blendSet
			if s in kickSet
				blendSet = setdiff(blendSet, s)
			end
		end
		println("[SBD] Kicking scenarios [$(length(kickSet))] $kickSet .")
		println("[SBD] Blending scneario [$(length(blendSet))] $blendSet .")

		if length(kickSet) >= ceil(exargs[:eps] * S)
			incumbentPool = setdiff(Set(collect(1:stoc.S)), kickSet)
			println("[II] Optimality reached :: enough scnearios are concluded to be definitely infeasible")
			stop = true
		end

		# Termination III :: Gap closing
		if gaptofixedLB - 0.0 < config.TOLERANCE
			incumbentPool = pickScenarioPool
			println("[III] Optimality reached :: closing gap to be 0.0")
			stop = true
		end

		# Terminate the algorithm
		if stop
			println("[POST] Printing columns into file. ")
			columns_printer(stoc)

			println("[POST] Considered scenarios $incumbentPool for POST-OPTIMIZATION")
			ubPool = incumbentPool
			ubISOCost = []
			for s in 1:S
				if s in ubPool
					push!(ubISOCost, isoCost[s])
				end
			end
			terminate, postProb = sbd_ubPool_solve(power, param, stoc, exargs, fixedLB, ubPool, ubISOCost, subprob_formulation, max(0.1,21600-tictoc))
			println("[POST] POST-OPTIMIZATION objective = ", getobjectivevalue(postProb.model), "[", getsolvetime(postProb.model), "s]")
			tictoc += getsolvetime(postProb.model)

			sbdSolution = get_primal_solution(postProb)
			attemptSBDDesign = get_design(postProb.model)
			if attemptSBDDesign.cost < incumbent
				sbdDesign = attemptSBDDesign
			else
				sbdDesign = incumbentDesign
			end

			totalcost, expandcost, hardencost = get_design_cost(sbdDesign, param)
			println(string("[POST] The total cost is ", totalcost, " = $expandcost + $hardencost"))
			println("[POST] Gap towards fixed LB [$(round.((totalcost-fixedLB)/totalcost*100,2))]")
			write_output_files(power, param, stoc, sbdDesign, exargs)
			print_design(sbdDesign, param)
			println("[EXIT] TicToc Now is $tictoc s")
			println("[POST] Successfully completed SBD heuristic.")
			flush(STDOUT)
			return postProb, sbdDesign
		end

		# =================================================================== #
		#         This section is for generating columns heuristic            #
		# =================================================================== #
		heuristicTime =
			@elapsed stoc, solved, tictoc =
				eval(parse(exargs[:CGHEURISTIC]))(power,
									param,
									stoc,
									exargs,
									masterDesign.cost,
									kickSet,
									blendSet,
									pickScenarioPool,
									neglectedScenarioPool,
									master_formulation,
									subprob_formulation,
									solved,
									allSubprobs,
									fixedLB = fixedLB,
									tictoc = tictoc,
									incumbent = incumbent,
									iterIdxer = iter,
									cuts = cuts)
		# ==================================================================== #
		# Reveal current columns infromation
		print_stoc_summary(stoc)
	end
end
