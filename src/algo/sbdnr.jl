function sbd_norisk(power::Dict, param::Dict, stoc::stocType, driver::Dict; selection=[], isoCost=[])

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

	infeas = [s for s in selection if incumbCol.feamap[s] == 0]

	while true
		iter += 1
		st = time()
		nextidx, infeas, feas = check_slackness(param, stoc, incumbCol, driver, selection=infeas, builtm=builtmodel)
		info("[SBD-NR] SLACKNESS $(time()-st) || A-T $(time()-ast)")

		cover = length(feas)/S*100
		info("[SBD-NR] incumbent cover rate $(length(feas)) [$(round.(cover,3))%]")

		if isempty(infeas) || iter > driver[:MAXITER]
			info("[SBD-NR: EXIT] Termination condition reached...")
			if !(driver[:SOLVER] == "Clp")
				if iter > 1
					solution = get_primal_solution(masterprob)
					write_output_files(power, param, stoc, solution, driver)
				else
					@show scenpool
					info("[SBD-NR: EXIT] Post optimization initiated on scenario set $(scenpool)")
					postprob = build_sp(param, stoc, driver, selection=scenpool, sbtype="tight")
					warmstart_heuristic(postprob, stoc, driver, selection=scenpool)
					config_solver(postprob, driver, timelimit=driver[:TIMELIMITII],focus="optimality")
					status = solve(masterprob.model)
					masterTime = getsolvetime(masterprob.model)
					tictoc += masterTime
					solution = get_primal_solution(masterprob)
				end
			end
			totalcost, expandcost, hardencost = get_design_cost(incumbCol, param)
			info("[SBD-NR:EXIT] The total cost is $totalcost = $expandcost + $hardencost")
			print_design(incumbCol, param)
			write_output_files(power,param,stoc,solution,driver)
			return masterprob, incumbCol
		end

		push!(scenpool, nextidx)	# For next master problem

		info("[SBD-NR] master pool $scenpool")
		st = time()
		masterprob = build_sp(param, stoc, driver, selection=scenpool, sbtype="tight")
		warmstart_heuristic(masterprob, param, stoc, driver, selection=scenpool)
		config_solver(masterprob, driver, timelimit=driver[:TIMELIMITII], focus="optimality")
		status = solve(masterprob.model)

		if status == :Infasible  # Assumption for climate problem
			print_iis_gurobi(masterprob.model)
			error("[SBD-NR] Master infeasible")
		end

		mBound = masterprob.model.objBound
		mObj = master.model.objVal
		mTime = time()-st
		incumbCol = get_design(masterprob, idx=iter, lb=mBound, time=mTime)
		info("[SBD-NR] ITER $(iter) || UB $(mObj) || LB $(mBound)")
		info("[SBD-NR] MASTER $(time()-st) || A-TIME $(time()-ast)")
	end
end
