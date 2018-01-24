function build_sp(param::Dict, stoc::stocType, selection, driver::Dict; sbtype="free")

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

function sbd_master_formulation(power::Dict, param::Dict, stoc::stocType, driver::Dict, overrideEps=driver[:eps], cuts=[], prevSolution=[]; kwargs...)

	options = Dict(kwargs)
	master = oneProblem()
	m = init_model_solver()	# setup empty model and associated solver

	# Attach the parameters
	master.param = param

	# Stoc always remains to be the full set
	S = stoc.S
	T = master.param[:T] = stoc.T
	B = master.param[:B]
	L = master.param[:L]

	master.name = power["name"]
	master.stage = 1
	master.T = stoc.T
	master.vars = Dict()

	# Master problem capture the scenarios that is assinged to it
	master.S = master.param[:S] = stoc.S

	# Capture the biggest solution pool
	poolLength = length(stoc.sbdColumns)

	master.vars[:pg] = @variable(m, 0<=pg[i=1:B, 1:T]<=param[:Pgbar][i], Int)
	master.vars[:h] = @variable(m, 0<=h[i=1:B, 1:T]<=param[:Hbar][i], Int)
	master.vars[:Y] = @variable(m, Y[k=1:poolLength], Bin)	# Design selector
	master.vars[:ACT] = @variable(m, ACT[1:S], Bin) 		# Indicates if a scenarios has been activated for feasible or not

	if isa(prevSolution, solnType)
		info("[SBD] Warm starting master problem with known solutions...")
		for i in 1:B
			for j in 1:T
				setvalue(master.vars[:pg][i,j], prevSolution.primal[:pg][i,j])
				setvalue(master.vars[:h][i,j], prevSolution.primal[:h][i,j])
			end
		end
		for i in 1:length(master.vars[:Y])
			if i <= length(prevSolution.primal[:Y])
				setvalue(master.vars[:Y][i], prevSolution.primal[:Y][i])
			else
				setvalue(master.vars[:Y][i], 0.0)
			end
		end
		for i in 1:length(master.vars[:ACT])
			setvalue(master.vars[:ACT][i], prevSolution.primal[:ACT][i])
		end
	end

	@constraint(m, [i=1:B], pg[i,1] >= param[:Pg0][i])
	@constraint(m, [i=1:B], h[i,1] >= param[:H0][i])
	@constraint(m, [i=1:B,t=2:T], pg[i,t-1] <= pg[i,t])
	@constraint(m, [i=1:B,t=2:T], h[i,t-1] <= h[i,t])

	@constraint(m, [i=1:B,t=1:T,k=1:poolLength; stoc.sbdColumns[k].active == true],
		pg[i,t]>=Y[k]*stoc.sbdColumns[k].pg[i,t])
	@constraint(m, [i=1:B,t=1:T,k=1:poolLength; stoc.sbdColumns[k].active == true],
		h[i,t]>=Y[k]*stoc.sbdColumns[k].h[i,t])

	@constraint(m, avoid, sum(Y[k] for k=1:poolLength if stoc.sbdColumns[k].active == false) == 0)

	@constraint(m, [s=1:S, k=1:poolLength; stoc.sbdColumns[k].active == true],
		ACT[s] >= Y[k]*stoc.sbdColumns[k].feamap[s])
	@constraint(m, [s=1:S],
		ACT[s] <= sum(Y[k]*stoc.sbdColumns[k].feamap[s] for k=1:poolLength))

	@constraint(m, risk,
		sum(ACT[s] for s=1:S) >= master.S * (1-overrideEps))

	if haskey(options, :incumbent)
		info("[ISO] Adding valid inequalities to the master problem")
		@constraint(m, cuts_expr[k=1:poolLength; stoc.sbdColumns[k].lb > options[:incumbent]],
			sum(ACT[s] for s in stoc.sbdColumns[k].source) <= length(stoc.sbdColumns[k].source)-1)
	end

	@objective(m, Min,
		sum(param[:Cg][i,1]*(pg[i,1]-param[:Pg0][i]) for i=1:B) + sum(param[:Ch][i,1]*(h[i,1]-param[:H0][i]) for i=1:B) + sum(param[:Cg][i,t]*(pg[i,t]-pg[i,t-1]) for i=1:B,t=2:T) + sum(param[:Ch][i,t]*(h[i,t]-h[i,t-1]) for i=1:B,t=2:T));

	master.model = m;

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
