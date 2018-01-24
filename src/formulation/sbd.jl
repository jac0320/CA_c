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

function check_feasible(param::Dict, stoc::stocType, soln, driver::Dict; selection=[], builtm=nothing)

	isempty(selection) ? selection = [1:param[:S];] : selection = selection
	feamap = zeros(Int, param[:S])
	for s in selection
		if builtm != nothing
			fix_xbar(builtm[s], param, soln)
			subprob = builtm[s]
		else
			subprob = sp_formulation(param, soln)
			attach_scenario(subprob, param, stoc, [s], driver[:MODEL], driver, sbtype="tight")
		end
		# reason_logical_vars(subprob, param, stoc, s)
		config_solver(subprob.model, driver, timelimit=20, focus="feasibility",threads=1)
		status = solve(subprob.model, suppress_warnings=true)
		status != :Infeasible ? feamap[s] = 1 : feamap[s] = 0
	end

	return feamap
end

function check_slackness(param::Dict, stoc::stocType, soln, driver::Dict; selection=[], builtm=nothing)

	feaCnt = 0
	feaDict = Dict()
	infeaDict = Dict()
	maxSlack = -Inf
	maxSlackIdx = 0
	coverage = 0.0

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	slackPool = zeros(length(selection))

	# Performing feasibility check on all scenarios
	for s in 1:length(selection)
		slackInd = pop!(selection)
		if builtm != nothing
			fix_xbar(builtm[slackInd], param, soln)
			slackProb = builtm[slackInd]
		else
			slackProb = sp_formulation(param, soln)
			attach_scenario(slackProb, param, stoc, [slackInd], driver[:MODEL], driver, sbtype="slackness")
		end
		# reason_logical_vars(slackProb, param, stoc, s)
		config_solver(slackProb.model, driver, timelimit=180, focus="optimality", threads=1)
		print(slackProb.model)
		error("STOP")
		status = solve(slackProb.model, suppress_warnings=true)
		status == :Infeasible && print_iis_gurobi(slackProb.model, driver)
		status == :Infeasible && error("Slackness problem should always be feasible...")
		push!(slackPool, getobjectivevalue(slackProb.model))
		abs(slackPool[end]) > 0.001 ? infeaDict[slackInd] = slackPool[end] : feaDict[slackInd] = slackPool[end]
		if slackPool[end] > maxSlack
			maxSlack = slackPool[end]
			maxSlackIdx = slackInd
		end
	end

	return maxSlackIdx, infeaDict, feaDict
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

	@assert size(xbar.pg) == (B,T)
	@assert size(xbar.h) == (B,T)

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
			expandVal = getupperbound(prob.vars[:pg][i,t])
			hardenVal = getupperbound(prob.vars[:h][i,t])
			ss = stoc.scenarios[s].data["SS"][i,t]
			maxharden = param[:ProM][i]*hardenVal + param[:Ele][i]
			ss - maxharden > 0 ? assDet[i,t] = 0 : assDet[i,t] = 1
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
