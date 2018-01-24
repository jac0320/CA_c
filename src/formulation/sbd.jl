function build_sp(param::Dict, stoc::stocType, driver::Dict; selection::Any=[], sbtype="free")

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	sp = sp_formulation(param)
	attach_scenario(sp, param, stoc, selection, driver[:MODEL], driver, sbtype=sbtype)

    return sp
end

function build_sp(prob::oneProblem, param::Dict, soln::Any)
	fix_xbar(prob, param, soln)
	return prob
end

function sp_formulation(param::Dict, soln=nothing)

	prob = oneProblem()
	post_adaptation_vars(prob, param)

	soln != nothing && fix_xbar(prob, param, soln)
	post_incremental_cons(prob, param)
	post_adaptation_obj(prob, param)

	return prob
end

function attach_scenario(prob::oneProblem, param::Dict, sto::stocType, selection, dispatch_model::Function, driver::Dict; sbtype="tight")

	isa(selection, Set) ? selection = collect(selection) : selection = selection
	isempty(selection) && error("setting up subproblem with 0 scenarios")

	post_logical_vars(prob, param, selection)
	if sbtype in ["tight", "slackness"]
		driver[:ALGO] != "evaluate" && post_logical_cons(prob, param, sto, selection)
	else
		post_logical_cons(prob, param, sto, selection)
	end
	post_risk_cons(prob, param, driver, selection, overrideEps=0.0)

    dispatch_model(prob, param, driver, selection, sbtype)

    return prob
end

function check_feasible(power::Dict, param::Dict, stoc::stocType, soln, driver::Dict, selection=[])

	infeaPool = []
	feaPool = []
	ub = Inf
	ub_coverage = 0.0

	isempty(selection) ? selection = [1:S;] : selection = selection

	if haskey(options, :builtModel)
		allSubprobs = options[:builtModel]
		@assert isa(allSubprobs, Array{oneProblem})
		@assert length(allSubprobs) == length(selection)
		builtModel = true
	else
		builtModel = false
	end

	for s in selection
		if builtModel
			@assert length(selection) == stoc.S
			fix_xbar(allSubprobs[s], param, soln)
			reason_logical_vars(allSubprobs[s], param, stoc, s)
			subprob = allSubprobs[s]
		else
			subprob = sp_formulation(param, soln)
			attach_scenario(subprob, param, stoc, [s], driver[:MODEL], driver, sbtype="tight")
			reason_logical_vars(subprob, param, stoc, s)
		end
		config_solver(subprob.model, driver, timelimit=20, focus="feasibility", threads=1)
		status = solve(subprob.model, suppress_warnings=true)
		status != :Infeasible ? push!(feaPool,s) : push!(infeaPool,s)
	end
	@assert length(feaPool) + length(infeaPool) == param[:S]

	ub_coverage = length(feaPool) / param[:S]
	ub = get_design_cost(soln, param)

	return ub, ub_coverage, infeaPool, feaPool
end

function check_slackness(power::Dict, param::Dict, stoc::stocType, soln, driver::Dict, selection=[]; sbtype="free")

	feaCnt = 0
	feaDict = Dict()
	infeaDict = Dict()
	maxSlack = -Inf
	maxSlackIdx = 0
	coverage = 0.0

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	slackPool = zeros(length(selection))

	if haskey(options, :builtModel)
		allSubprobs = options[:builtModel]
		@assert isa(allSubprobs, Array{oneProblem})
		@assert length(allSubprobs) >= length(selection) # All problems must be presented here
		builtModel = true
	else
		builtModel = false
	end

	# Performing feasibility check on all scenarios
	for s in 1:length(selection)
		slackInd = pop!(selection)
		if builtModel
			fix_xbar(allSubprobs[slackInd], param, soln)
			slackProb = allSubprobs[slackInd]
		else
			slackProb = sp_formulation(param, soln)
			attach_scenario(slackProb, param, stoc, [slackInd], driver[:MODEL], driver, sbtype="slackness")
			reason_logical_vars(slackProb, param, stoc, s)
		end
		reason_logical_vars(slackProb, param, stoc, s)
		config_solver(slackProb.model, driver, timelimit=180, focus="optimality", threads=1)
		status = solve(slackProb.model, suppress_warnings=true)
		status == :Infeasible && error("Slackness problem should always be feasible...")
		push!(slackPool, getobjectivevalue(slackProb.model))
		abs(slackPool[end]) > 0.001 ? infeaDict[slackInd] = slackPool[end] : feaDict[slackInd] = slackPool[end]
		if slackPool[end] > maxSlack
			maxSlack = slackPool[end]
			maxSlackIdx = slackInd
		end
	end

	return maxSlackIdx, maxSlack, infeaDict, feaDict
end

function sbd_master_formulation(power::Dict, param::Dict, stoc::stocType, driver::Dict, overrideEps=driver[:eps], cuts=[], prevSol=[]; incumbent=Inf)

	master = oneProblem()

	post_adaptation_vars(master, param)
	post_master_logical_vars(master, param, stoc.sbdColumns)
	post_incremental_cons(master, param)

	post_union_operator(master, param, stoc.sbdColumns)
	post_column_onoff(master, stoc.sbdColumns)
	post_scenario_activation(master, param, stoc.sbdColumns)
	post_risk_con(master, param[:S], driver, :ACT)
	post_master_cuts(master, stoc.sbdColumns, incumbent)

	if isa(prevSol, solnType)
		prevYval = [prevSol.primal[:Y],zeros(P-length(prevSol.primal[:Y]));]
		enforce_startval(master, :pg, sval=prevSol.primal[:pg])
		enforce_startval(master, :h, sval=prevSol.primal[:h])
		enforce_startval(master, :Y, sval=prevYval)
		enforce_startval(master, :ACT, sval=prevSol.primal[:ACT])
	end

	post_adaptation_obj(master, param)

	return master
end

function fix_xbar(prob::oneProblem, param::Dict, xbar::solnType)

	B, T = param[:B], param[:T]

	@assert size(x_bar.primal[:pg]) == (B,T)
	@assert size(x_bar.primal[:h]) == (B,T)

	for i=1:B
		for j=1:T
			setupperbound(prob.vars[:pg][i,j], xbar.primal[:pg][i,j])
			setlowerbound(prob.vars[:pg][i,j], xbar.primal[:pg][i,j])
			setupperbound(prob.vars[:h][i,j], xbar.primal[:h][i,j])
			setlowerbound(prob.vars[:h][i,j], xbar.primal[:h][i,j])
		end
	end

	return
end

function fix_xbar(prob::oneProblem, param::Dict, xbar::designType)

	B, T = param[:B], param[:T]

	@assert size(x_bar.pg) == (B,T)
	@assert size(x_bar.h) == (B,T)

	for i=1:B
		for j=1:T
			setupperbound(prob.vars[:pg][i,j], xbar.pg[i,j])
			setlowerbound(prob.vars[:pg][i,j], xbar.pg[i,j])
			setupperbound(prob.vars[:h][i,j], xbar.h[i,j])
			setlowerbound(prob.vars[:h][i,j], xbar.h[i,j])
		end
	end

	return
end

function reason_logical_vars(prob::oneProblem, param::Dict, stoc::stocType, s::Int)

	B, T = param[:B], param[:T]
	assDet = ones(Int, B, T)

	for i in 1:B
		for t in 1:T
			expandVal = getupperbound(subprob.vars[:pg][i,t])
			hardenVal = getupperbound(subprob.vars[:h][i,t])
			ss = stoc.scenarios[s].data["SS"][i,t]
			maxharden = param[:ProM][i]*hardenVal + param[:Ele][i]
			ss - maxharden > 0 ? subprob.param[:assDet][i,t,1] = 0 : subprob.param[:assDet][i,t,1] = 1
			@assert	haskey(param, :aslDet)
			setupperbound(prob.vars[:a][i,t,1], param[:aslDet][i,t,s] * assDet[i,t])
			setlowerbound(prob.vars[:a][i,t,1], param[:aslDet][i,t,s] * assDet[i,t])
			setupperbound(prob.vars[:ap][i,t,1], param[:aslDet][i,t,s] * assDet[i,t] * expandVal)
			setlowerbound(prob.vars[:ap][i,t,1], param[:aslDet][i,t,s] * assDet[i,t] * expandVal)
			setupperbound(prob.vars[:ah][i,t,1], assDet[i,t] * hardenVal)
			setlowerbound(prob.vars[:ah][i,t,1], assDet[i,t] * hardenVal)
			setupperbound(prob.vars[:ass][i,t,1], assDet[i,t])
			setlowerbound(prob.vars[:ass][i,t,1], assDet[i,t])
		end
	end

	return
end
