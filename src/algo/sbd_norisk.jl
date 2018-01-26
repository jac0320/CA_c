function sbd_norisk(power::Dict, param::Dict, stoc::stocType, exargs::Dict,
					subprob_formulation::Function,
					selection=[], isoCost=[]; kwargs...)

	options = Dict(kwargs)

	isempty(selection) ? S = stoc.S : S = length(selection)
	(isempty(selection)) && (selection = [1:S;])
	(haskey(options, :isoCostPool)) && (isoCostPool = options[:isoCostPool])

	iter = 0
	tictoc = 0.0
	masterPool = Int[]
	feaPool = Set()
	infeaPool = Set()

	# if isempty(isoCost)
	# 	println("[SBD-NORISK] ISO stage ...")
	# 	isoTime, stoc, isoUnionCost, isoCost, earlyExit =
	# 		isolate_stage(power, param, stoc, exargs, subprob_formulation, selection)
	# 	tictoc += isoTime
	# 	println("[TICTOC] After ISO stage = $tictoc s")
	# end
    #

	if isempty(isoCost)
		initScenIdx = 1
		incumbentDesign = isolate_solve_one_scenario(power, param, stoc, initScenIdx, exargs, subprob_formulation, "free")
		incumbentDesign.feamap = zeros(Int, S)
		push!(masterPool, initScenIdx)
	else
		initScenIdx = findfirst(isoCost, maximum(isoCost))
		incumbentDesign = stoc.sbdColumns[initScenIdx]
		push!(masterPool, initScenIdx)
		for s in selection
			(incumbentDesign.feamap[s]==0) ? push!(infeaPool, s) : push!(feaPool, s)
		end
	end

	allSubprobs = Array{oneProblem}(stoc.S)
	println("[SBD-NORISK] creating slackness problems...")
	for ss = 1:stoc.S
		allSubprobs[ss] = oneProblem()
		allSubprobs[ss] = sbd_base_formulation(power, param, stoc)
		allSubprobs[ss] = attach_scenario(allSubprobs[ss], stoc, [ss], exargs[:MODEL], 0.0, exargs, subprobType="slackness")
	end

	println("[SBD-NORISK] starting main loop...")
	stop = false

	while true
		iter += 1
		println("[SBD-NORISK] entering iteration $iter...")
		csTime = @elapsed maxSlackIdx, maxSlack, infeaDict, feaDict =
			check_slackness(power, param, stoc, incumbentDesign, exargs, infeaPool, builtModel=allSubprobs)
		tictoc += csTime
		println("[TICTOC] after checking slackness in iteration $iter TIME = $tictoc")

		newCover = length(feaDict)/(stoc.S)*100
		println("[SBD-NORISK] incumbent cover rate $(length(feaDict)) [$(round.(newCover,3))%]")

		if isempty(infeaDict) || iter == config.MAXITER || stop
			println("[SBD-NORISK: EXIT] No infeasible scenarios located. Termination condition reached...")
			println("[TICTOC] final time mark is $tictoc")
			if !(config.SOLVER == "Clp")
				if iter > 1
					solution = get_primal_solution(masterProb)
					write_output_files(power,param,stoc,solution,exargs)
				else
					println("[SBD-NORISK: EXIT] ??? Post optimization initiated on scenario set $(masterPool)")
					masterProb = subprob_formulation(power, param, stoc, masterPool, exargs)
					warmstart_heuristic(masterProb, power, param, stoc, exargs, selection=masterPool)
					solver_config(masterProb.model, timelimit=config.TIMELIMITII, mipgap=config.OPTGAP, showlog=1, focus="optimality", threads=16)
					status = solve(masterProb.model)
					masterTime = getsolvetime(masterProb.model)
					tictoc += masterTime
					println("[SBD-NORISK] 1 iteration master solved. NOW TIME = $tictoc")
					println("[SBD-NORSIK] Incumbent cost is $(get_design_cost(get_design(masterProb.model), param))")
					solution = get_primal_solution(masterProb)
					write_output_files(power, param,stoc,solution,exargs)
				end
			end
			totalcost, expandcost, hardencost = get_design_cost(incumbentDesign, param)
			println(string("[SBD-NORISK:EXIT] The total cost is ", totalcost, " = $expandcost + $hardencost"))
			print_design(incumbentDesign, param)
			return masterProb, incumbentDesign
		end

		push!(masterPool, maxSlackIdx)

		println("[SBD-NORISK] iteration $(iter) current master scenario pool $masterPool")
		masterProb = subprob_formulation(power, param, stoc, masterPool, exargs, subprobType="tight")
		warmstart_heuristic(masterProb, power, param, stoc, exargs, selection=masterPool)
		solver_config(masterProb.model, timelimit=config.TIMELIMITII, mipgap=config.OPTGAP, focus="optimality", showlog=config.SHOWLOG, threads=16)
		status = solve(masterProb.model)
		masterTime = getsolvetime(masterProb.model)
		masterBound = getobjbound(masterProb.model)
		println("[SBD-NORISK] This iteration's master solving time = $masterTime [LB=$(masterBound)]")
		tictoc += masterTime
		println("[TICTOC] after master problem in iteration $iter TIME = $tictoc")

		if status == :Infasible  # Assumption for climate problem
			print_iis_gurobi(masterProb.model)
			error("ERROR|sbd.jl|sbd_norisk()|Joint problem infeasible. Check the iis report on the problem.")
		end

		incumbentDesign = get_design(masterProb.model)
		incumbentDesign.cost = getobjectivevalue(masterProb.model)
		incumbentDesign.time = getsolvetime(masterProb.model)

		infeaPool = Set()
		[push!(infeaPool, s) for s in keys(infeaDict)]
		println("[SBD-NORISK] incumbent cost = $(incumbentDesign.cost)")
	end
end
