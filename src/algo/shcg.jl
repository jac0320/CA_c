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
	builtmodel = options[:builtmodel]

	jobs = shcg_generate_jobs(param, solved, dim=iter)
	info("[ENU] Job Count = $(length(jobs))")

	shcg_process_jobs(param, stoc, driver, builtmodel, incumb)
	# if driver[:PARALLEL]
	# 	tSolve = @elapsed cols = pmap((a1,a2,a3,a4)->shcg_solve_sp(a1,a2,a3,a4),
	# 									[param for j in jobs],
	# 									[stoc for j in jobs],
	# 									[driver for j in jobs],
	# 									[j for j in jobs])
	# 	info("[ENU] Solve time $(tsolve)")
    #
	# 	groups = shcg_regroup_cols(cols, driver)
	# 	tcheck = @elapsed groups = pmap((a1,a2,a3,a4)->shcg_check_fea(a1,a2,a3,a4),
	# 									[param for g in groups],
	# 									[stoc for g in groups],
	# 									[driver for g in groups],
	# 									[g for g in groups])
	# 	info("[ENU] Feasibility check time $(tcheck)")
    #
	# 	tstore = @elapsed [shcg_store_col(stoc.sbdColumns, g, incumb=incumb) for g in groups]
	# 	info("[ENU] Storaging time $(tstore)")
	# else
	# 	tsolve = @elapsed cols = [shcg_solve_sp(param, stoc, driver, j) for j in jobs]
	# 	info("[ENU] Solve time $(tsolve)")
    #
	# 	tcheck = @elapsed cols = [iso_check_block(c, param, stoc, driver, builtmodel=builtmodel) for c in cols]
	# 	info("[ENU] Feasibility check time $(tcheck)")
    #
	# 	tstore = @elapsed [shcg_store_col(stoc.sbdColumns, c, incumb=incumb) for c in cols]
	# 	info("[ENU] Storaging time $(tstore)")
	# end

	return
end

function shcg_enu(driver::Dict; findparameter=false)

	if findparameter
		if driver[:CGMAX] > 0
			return Dict(:noimprovestop=>driver[:CGMAX])
		else
			return Dict(:noimprovestop=>6)
		end
	end

	return Dict
end

function shcg_heu(param::Dict, stoc::stocType, driver::Dict; kwargs...)

	options = Dict(kwargs)

	iter = options[:iter]
	incumb = options[:incumb]
	solved = options[:solved]
	isoCost = options[:isoCost]
	builtmodel = options[:builtmodel]

	S = param[:S]

	ranklist = sort(collect(zip(isoCost,[1:S;])))

	pivotrank = mod(iter, param[:S]) + 1
	dim = convert(Int, floor(iter/param[:S]))
	pivot = ranklist[pivotrank][2]

	# Pivot with selected job
	jobs = shcg_generate_jobs(param, solved, pivot, dim=dim)

	info("Running CG heuristic :: lineup_partial_heuristic")
	incumblb = extras[:incumblb]

	shcg_process_jobs(param, stoc, driver, builtmodel, incumb)

	return
end

function shcg_nr(param::Dict, stoc::stocType, driver::Dict; selection=[], isoCost=[])

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
