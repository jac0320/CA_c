function shcg_solve_master(param::Dict, stoc::stocType, driver::Dict; prevsol=nothing, incumb=Inf, builtmodel=nothing)

	master = sbd_master_formulation(param, stoc, driver, prevsol=prevsol, incumbent=incumb)
	config_solver(master.model, driver, timelimit=driver[:TIMELIMITII])
	status = solve(master.model, suppress_warnings=true)

	if status == :Infeasible
		print_iis_gurobi(master.model)
		error("ERROR|alg_sbd.jl|shcg()|master proble infeasible, check formulation");
	else
		mastercol = get_design(master)
		mastercol.feamap = check_feasible(param, stoc, mastercol, driver, builtmodel=builtmodel)
		mastercol.coverage = sum(mastercol.feamap) / param[:S]
		sum(mastercol.feamap) < (1-param[:eps]) * param[:S] && warn("POTENTIAL BUG :: master solution problematic !!")
		shcg_store_col(stoc.sbdColumns, mastercol, incumb=Inf)
	end

	return get_primal_solution(master)
end

function shcg_solve_sp(param::Dict, stoc::stocType, driver::Dict, selection::Vector=[])

	onesp = build_sp(param, stoc, driver, selection=selection, sbtype="tight")
	config_solver(onesp.model, timelimit=driver[:TIMELIMITIII], focus="optimality", threads=driver[:WORKERTHREADS])
	status = solve(onesp.model, suppress_warnings=true)

	if status == :Infeasible
		print_iis_gurobi(onesp.model)
		error("ERROR|alg_sbd_jl|shcg|Subprob infeasible, check formulation")
	end

	col = get_design(onesp)

	return col
end

function shcg_solve_sp(param::Dict, stoc::stocType, driver::Dict, selection::Any, isoCost::Array)

	col = shcgnr(param, stoc, driver, selection=selection, isoCost=isoCost)
	col.feamap = check_feasible(param, stoc, d, driver)
	col.coverage = sum(col.feamap) / param[:S]

	return col
end

function shcg_check_fea(param::Dict, stoc::stocType, driver::Dict, colGroup::Array{designType})

	S = param[:S]

	ms = [build_sp(param, stoc, driver, selection=[s], sbtype="tight") for s in 1:S]

	for d in colGroup
		d.feamap = check_feasible(param, stoc, d, driver, builtmodel=ms)
		d.coverage = sum(d.feamap)/S
	end

	return colGroup
end

function search_incumb_col(pool::Vector{designType}, param::Dict)

	P = length(pool)

	incumbub = Inf

	for i in 1:P
		if pool[i].cost < incumbub && pool[i].coverage >= 1 - param[:eps]
			incumbub = pool[i].cost
			incumbcol = pool[i]
		end
	end

	return incumbcol
end

# TODO: add more
function shcg_generate_jobs(param::Dict, solved::Dict; dim::Int=1)

	S = param[:S]
	rotations = [i for i in 1:S]

	jobs = []

	for i in subsets(rotations, dim)
		Set(i) in solved || push!(jobs, i)
		Set(i) in solved || push!(solved, Set(i))
	end

	return jobs
end

function partial_generate_pool(param::Dict, solved::Dict, pivot::Int; dim::Int=1)

	S = param[:S]
	rotations = [i for i in 1:S if i != pivot]

	jobs = []

	for i in subsets(rotations, dim)
		j = [pivot;i]
		@assert length(Set(j)) == length(j)
		Set(j) in solved || push!(jobs, j)
		Set(j) in solved || push!(solved, Set(j))
	end

	return jobs
end

function shcg_store_col(pool::Array{designType}, design::designType; incumb=Inf)

	P = length(pool)
	same = false
	design.cost > incumb && return

	for i in 1:P
		if (design.pg == pool[i].pg) && (design.h == pool[i].h)
			same = true
			break
		end
	end

	design.k = length(pool) + 1
	push!(designPool, design)

    return
end

function shcg_store_col(pool::Array{designType}, design::Vector{designType}; incumb=Inf)

	for d in design
		shcg_store_col(pool, d, incumb=incumb)
	end

	return
end

function shcg_regroup_cols(pool::Array{designType}, driver::Dict)

	P = length(pool)
	W = driver[:WORKERS]

	rCnt = mod(P, W)
	groups =[]
	for i in partition(pool, W)
		push!(groups, i)
	end

	last_group = pool[end-rCnt+1:end]
	isempty(last_group) || push!(groups, last_group)

	return groups
end

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

		dPool = pmap((a1,a2,a3,a4,a5,a6)->iso_check_block(a1,a2,a3,a4,a5),
							[param for s in dPool],
							[stoc for s in dPool],
							[driver for d in dPool],
							[d for d in dPool],
							[S for d in dPool])

	else  # sequential implementation
		sp_models = [build_sp(param, stoc, driver, selection=[s], sbtype="tight") for s in selection]
		dPool = [iso_solve_block(sp_models[s], driver, source=[s]) for s in selection]
		for d in dPool
			iso_check_block(d, param, stoc, driver, builtmodel = sp_models)
		end
	end

	for s in selection
		collect_design(stoc.sbdColumns, dPool[s])
	end

	if !skipsafety
		ud = union_design(dPool, param, source=[1:S;])
		collect_design(stoc.sbdColumns, ud)
	end

	return
end

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

function iso_check_block(d::designType, param::Dict, stoc::stocType, driver::Dict; selection=[], builtmodel=nothing)

	d.feamap = check_feasible(param, stoc, d, driver, builtmodel=builtmodel)
	d.coverage = sum(d.feamap) / S

	return d
end


function columns_manager(stoc::stocType, kickScen, incumbent::Float64, driver::Dict; kwargs...)

	allScen = Set()
	for s in 1:stoc.S
		push!(allScen, s)
	end

	if isa(kickScen, Array)
		kickScenSet = Set(kickScen)
	elseif isa(kickScen, Set)
		kickScenSet = kickScen
	end
	activeScen = setdiff(allScen, kickScenSet)

	eliminateScen = Set()
	scenEliminater = Dict(i=>Set() for i=1:stoc.S)

	S = stoc.S - length(kickScen)
	requireS = ceil(stoc.S * (1-driver[:eps]))

	# First Summarize All Generated Columns
	for i in 1:length(stoc.sbdColumns)
		if incumbent + config.TOLERANCE < stoc.sbdColumns[i].lb && (length(S) < stoc.S) && length(stoc.sbdColumns[i].source) <= 2
			if stoc.sbdColumns[i].active
				println("[COLUMN MANAGER] Found eliminating column [scenario set $(stoc.sbdColumns[i].source)][Cost $(stoc.sbdColumns[i].cost)]\n")
			end
			stoc.sbdColumns[i].active = false
			# Inference phase
			for s in stoc.sbdColumns[i].source
				avoidSet = setdiff(Set(stoc.sbdColumns[i].source),s)
				if !isempty(avoidSet)
					union(scenEliminater[s], unique(avoidSet))
				end
			end
		end
	end

	for i in activeScen
		combineAvoid = collect(unique(scenEliminater[i]))
		if length(combineAvoid) > stoc.S - requireS
			push!(eliminateScen, i)
		end
	end

	if !isempty(eliminateScen)
		info("[COLUMN MANAGER] Helped eliminate $(length(eliminateScen))")
	end

	return collect(eliminateScen)
end

function columns_cleaner(stoc::stocType, incumbent::Float64; kwargs...)

	for i in 1:length(stoc.sbdColumns)
		if stoc.sbdColumns[i].lb >= incumbent + config.TOLERANCE && stoc.sbdColumns[i].active
			info("Deactivating column $(i) with cost $(round.(stoc.sbdColumns[i].cost,2)). [$(round.(incumbent,2))]\n")
			stoc.sbdColumns[i].active = false
		end
	end

	return 0
end
