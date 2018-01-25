function isolate_stage(param::Dict, stoc::stocType, driver::Dict; selection=[], sbtype="tight", skipsafety=false)

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	S = length(selection)

	# =========================== Isolation Stage ===========================#
	dPool = Array{designType}(S)

	if driver[:PARALLEL] # parallel implementation
		dPool = pmap((a1,a2,a3,a4,a5,a6)->iso_solve_block(a1,a2,a3,a4,a5,a6),
							[param for s in selection],
							[stoc for s in selection],
							[s for s in selection],
							[driver for s in selection],
							[isoprob_formulation for s in selection],
							[sbtype for s in selection])

		dPool = pmap((a1,a2,a3,a4,a5,a6)->fea_check_block(a1,a2,a3,a4,a5),
							[param for s in dPool],
							[stoc for s in dPool],
							[driver for d in dPool],
							[d for d in dPool],
							[S for d in dPool])

	else  # sequential implementation
		sp_models = [build_sp(param, stoc, driver, selection=[s], sbtype="tight") for s in selection]
		dPool = [iso_solve_block(sp_models[s], driver, source=[s]) for s in selection]
		for d in dPool
			fea_check_block(d, param, stoc, driver, builtm = sp_models)
		end
	end

	# info("[ISO] Scen $s -> cost = [$(round.(dPool[s].cost,2))][LB=$(round.(dPool[s].lb,2))][TIME=$(round.(dPool[s].time,2))][Cover $(round.(dPool[s].coverage,2))][CHECKFEATime $(round.(oneCheckFeaTime,2))]")
	for s in selection
		collect_design(stoc.sbdColumns, dPool[s])
	end

	# Early exit rule :: if minimum cost isolated scenario fits the risk constraint, then find feasible lower bound.
	if !skipsafety
		ud = union_design(dPool, param, source=[1:S;])
		collect_design(stoc.sbdColumns, ud)
	end

	return
end


# =========================================== #


function iso_solve_block(param::Dict, stoc::stocType, s::Int, driver::Dict,sbtype="tight")

	p = build_sp(param, stoc, driver, selection=[s], sbtype=sbtype)
	# warmstart_heuristic(p, stoc, driver, selection=[s])
	config_solver(p.model, driver, timelimit=driver[:TIMELIMITIII], focus="optimality", threads=driver[:WORKERTHREADS])
	status = solve(oneIsoProb.model)

	if status == :Infeasible
		info("Creating infea column on scenario $s")
		d = infea_design(s, param)
	else
		d = get_design(oneIsoProb.model)
		d.k = s
		d.source = [s]
		d.time = getsolvetime(oneIsoProb.model)
		d.lb = oneIsoProb.model.objBound
		d.active = true
	end

	return d
end

function iso_solve_block(prob::oneProblem, driver::Dict; source=[])

	config_solver(prob.model, driver, timelimit=driver[:TIMELIMITIII], focus="optimality", threads=driver[:WORKERTHREADS])
	status = solve(prob.model, suppress_warnings=true)
	if status == :Infeasible
		info("Creating infea column on scenario $s")
		return infea_design(s, param)
	else
		d = get_design(prob.model)
		d.source = source
		d.time = getsolvetime(prob.model)
		d.lb = prob.model.objBound
		d.active = true
	end

	return d
end

function fea_check_block(d::designType, param::Dict, stoc::stocType, driver::Dict; selection=[], builtm=nothing)
	d.feamap = check_feasible(param, stoc, d, driver, builtm=builtm)
	return
end
