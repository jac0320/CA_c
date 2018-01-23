"""
	Input ::
		power::Dict, param::Dict, stocData::stochType, soln::designType/solnType
	Usage ::
		This function builds a base problem to get ready for diffferent needs in the sample-based algorithms.
	Output ::
		prob::oneProblem
"""
function sbd_base_formulation(power::Dict, param::Dict, stocData::stocType, soln=nothing, warmstart=nothing, builtModel=nothing)

	S = 0					#A clean master problem shouldn't carry any scenario
	T = param[:T]
	B = param[:B]

	if builtModel != nothing
		if isa(soln,solnType) == true
			for i=1:B
				for j=1:T
					setupperbound(builtModel.vars[:pg][i,j], soln.primal[:pg][i,j])
					setlowerbound(builtModel.vars[:pg][i,j], soln.primal[:pg][i,j])
					setupperbound(builtModel.vars[:h][i,j], soln.primal[:h][i,j])
					setlowerbound(builtModel.vars[:h][i,j], soln.primal[:h][i,j])
				end
			end
		elseif isa(soln,designType) == true
			for i=1:B
				for j=1:T
					setupperbound(builtModel.vars[:pg][i,j], soln.pg[i,j])
					setlowerbound(builtModel.vars[:pg][i,j], soln.pg[i,j])
					setupperbound(builtModel.vars[:h][i,j], soln.h[i,j])
					setlowerbound(builtModel.vars[:h][i,j], soln.h[i,j])
				end
			end
		else
			error("ERROR|climate_model.jl|sbd_base_formulation|Unkown solution type");
		end
		return builtModel
	else
		prob = oneProblem()
		prob.model = init_model_solver()	# setup empty model and associated solver
		prob.param = deepcopy(param);
		prob.vars = Dict()

		prob.name = power["name"]
		prob.T = stocData.T
		prob.S = 0			# Basic master problem is clean

		# add geneartion decision variables
		if soln == nothing
			prob.vars[:pg] = @variable(prob.model, 0<=pg[i=1:B, 1:T]<=param[:Pgbar][i], Int);
			prob.vars[:h] = @variable(prob.model, 0<=h[i=1:B, 1:T]<=param[:Hbar][i], Int);
		elseif isa(soln,solnType) == true
			prob.vars[:pg] = @variable(prob.model, 0<=pg[i=1:B, 1:T]<=param[:Pgbar][i], Int);
			prob.vars[:h] = @variable(prob.model, 0<=h[i=1:B, 1:T]<=param[:Hbar][i], Int);
			for i=1:B
				for j=1:T
					setupperbound(prob.vars[:pg][i,j], soln.primal[:pg][i,j])
					setlowerbound(prob.vars[:pg][i,j], soln.primal[:pg][i,j])
					setupperbound(prob.vars[:h][i,j], soln.primal[:h][i,j])
					setlowerbound(prob.vars[:h][i,j], soln.primal[:h][i,j])
				end
			end
		elseif isa(soln,designType) == true
			prob.vars[:pg] = @variable(prob.model, 0<=pg[i=1:B, 1:T]<=param[:Pgbar][i], Int);
			prob.vars[:h] = @variable(prob.model, 0<=h[i=1:B, 1:T]<=param[:Hbar][i], Int);
			for i=1:B
				for j=1:T
					setupperbound(prob.vars[:pg][i,j], soln.pg[i,j])
					setlowerbound(prob.vars[:pg][i,j], soln.pg[i,j])
					setupperbound(prob.vars[:h][i,j], soln.h[i,j])
					setlowerbound(prob.vars[:h][i,j], soln.h[i,j])
				end
			end
		else
			error("ERROR|climate_model.jl|sbd_base_formulation|Unkown solution type");
		end

		# (First Stage) Incremental Design Geneartion || Hardning
		@constraint(prob.model, [i=1:B], pg[i,1] >= param[:Pg0][i]);
		@constraint(prob.model, [i=1:B], h[i,1] >= param[:H0][i]);
		@constraint(prob.model, [i=1:B,t=2:T], pg[i,t-1] <= pg[i,t]);
		@constraint(prob.model, [i=1:B,t=2:T], h[i,t-1] <= h[i,t]);

		@objective(prob.model, Min, sum(param[:Cg][i,1]*(pg[i,1]-param[:Pg0][i]) for i=1:B)
							+ sum(param[:Ch][i,1]*(h[i,1]-param[:H0][i]) for i=1:B)
							+ sum(param[:Cg][i,t]*(pg[i,t]-pg[i,t-1]) for i=1:B,t=2:T)
							+ sum(param[:Ch][i,t]*(h[i,t]-h[i,t-1]) for i=1:B,t=2:T));
	end

	return prob
end

"""
	Input ::
		power::Dict, stoc::stocType, selec::Array/Set, characteristic::Function, overrideEps::Float64, exargs::Dict
	Usage ::
		This function takes a half-established base problem from sbd_base_formulation and expand the formulation by adding extra scenarios. The number of scenarios is flexible in this subroutine. It is widely used in sbd algorithms for constructing subproblems.
	Output ::
		prob::oneProblem
"""
function attach_scenario(prob::oneProblem, stoc::stocType, select, characteristic::Function, overrideEps, exargs::Dict; kwargs...)

	options = Dict(kwargs)

	if haskey(options, :subprobType)
	 	subprobType = options[:subprobType]
	else
		subprobType = "free"
	end

    # Measure the number of scenarios to be attahced
    S = length(select)

	# If no scenario is attached, simply return the original input
	if S == 0
		return prob
	end

	# Update the parameters in different structures
    prob.S = S			# Parameter updates
	prob.param[:S] = S 	# Parameter updates
    T = prob.T
    B = prob.param[:B]

    # Fetch some existing variables
	@assert haskey(prob.vars, :pg)
	@assert haskey(prob.vars, :h)

    # Introduce the necessary variables
	prob.vars[:f] = @variable(prob.model, f[1:T, 1:S], Bin)
	prob.vars[:a] = @variable(prob.model, a[1:B, 1:T, 1:S], Bin)
	prob.vars[:ass] = @variable(prob.model, ass[1:B, 1:T, 1:S], Bin)
	prob.vars[:ap] = @variable(prob.model, ap[1:B, 1:T, 1:S], Int)
	prob.vars[:ah] = @variable(prob.model, ah[1:B, 1:T, 1:S], Int)

	# If scenario selection is a type set, need to conver it into array
	if isa(select, Set) == true
		tempArrray = []
		for i in 1:length(select)
			push!(tempArray, pop!(select))
		end
		# Replace the set input with a translated array
		select = tempArray
	end

	if subprobType == "tight" || subprobType == "slackness"

		# prob.param[:assDet] = ones(Int, B, T, S)
		# for i in 1:B
		# 	for t in 1:T
		# 		for s in 1:S
		# 			hardenVal = getupperbound(prob.vars[:h][i,t])
		# 			expandVal = getupperbound(prob.vars[:pg][i,t])
		# 			if stoc.scenarios[select[s]].data["SS"][i,t] - (prob.param[:ProM][i]*hardenVal + prob.param[:Ele][i]) >= 0
		# 				prob.param[:assDet][i,t,s] = 0
		# 			end
		# 			@constraint(prob.model,
		# 				prob.vars[:a][i,t,s] == prob.param[:aslDet][i,t,select[s]] * prob.param[:assDet][i,t,s])
		# 			@constraint(prob.model,
		# 				prob.vars[:ap][i,t,s] == prob.vars[:a][i,t,s] * expandVal)
		# 			@constraint(prob.model,
		# 				prob.vars[:ah][i,t,s] == prob.param[:assDet][i,t,s] * hardenVal)
		# 			@constraint(prob.model,
		# 				prob.vars[:ass][i,t,s] == prob.param[:assDet][i,t,s])
		# 		end
		# 	end
		# end

		# Sea Level | Storm Surge Feasible Indicator
		if exargs[:ALGO] != "evaluate"
			@constraint(prob.model, [i=1:B, t=1:T, s=1:S],
				stoc.scenarios[select[s]].data["SS"][i,t] * (2*prob.vars[:ass][i,t,s]-1) - 2*prob.param[:Ele][i]*prob.vars[:ass][i,t,s] - 2*prob.param[:ProM][i]*prob.vars[:ah][i,t,s] + prob.param[:Ele][i] + prob.param[:ProM][i]*prob.vars[:h][i,t] <= 0)

			for i in 1:B
				for t in 1:T
					for s in 1:S
						@assert	haskey(prob.param, :aslDet)
						@constraint(prob.model, prob.vars[:a][i,t,s] == prob.param[:aslDet][i,t,select[s]] * prob.vars[:ass][i,t,s])
						mccormick(prob.model, prob.vars[:ap][i,t,s], prob.vars[:a][i,t,s], prob.vars[:pg][i,t], 0, 1, prob.param[:Pg0][i], prob.param[:Pgbar][i])
						mccormick(prob.model, prob.vars[:ah][i,t,s], prob.vars[:ass][i,t,s], prob.vars[:h][i,t], 0, 1, prob.param[:H0][i], prob.param[:Hbar][i])
					end
				end
			end
		end

	else

		# Sea Level | Storm Surge Feasible Indigcator
		@constraint(prob.model, [i=1:B, t=1:T, s=1:S],
			stoc.scenarios[select[s]].data["SS"][i,t] * (2*prob.vars[:ass][i,t,s]-1) - 2*prob.param[:Ele][i]*prob.vars[:ass][i,t,s] - 2*prob.param[:ProM][i]*prob.vars[:ah][i,t,s] + prob.param[:Ele][i] + prob.param[:ProM][i]*prob.vars[:h][i,t] <= 0)

	    for i in 1:B
	        for t in 1:T
	            for s in 1:S
					@assert	haskey(prob.param, :aslDet)
					@constraint(prob.model, prob.vars[:a][i,t,s] == prob.param[:aslDet][i,t,select[s]] * prob.vars[:ass][i,t,s])
	                mccormick(prob.model, prob.vars[:ap][i,t,s], prob.vars[:a][i,t,s], prob.vars[:pg][i,t], 0, 1, prob.param[:Pg0][i], prob.param[:Pgbar][i])
					mccormick(prob.model, prob.vars[:ah][i,t,s], prob.vars[:ass][i,t,s], prob.vars[:h][i,t], 0, 1, prob.param[:H0][i], prob.param[:Hbar][i])
	            end
	        end
	    end

	end

    characteristic(prob, exargs, select, subprobType=subprobType)
    formulation_features(prob, stoc, exargs, overrideEps, select, subprobType=subprobType);

    return prob
end

function check_feasible(power::Dict, param::Dict, stoc::stocType, soln, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)

	if haskey(options, :warmStarter)
		@assert isa(options[:warmStarter], solnType)
		warmstart = true
	else
		warmstart = false
	end

	feaCnt = 0
	infeaPool = []
	feaPool = []
	ub = Inf
	ub_coverage = 0.0

	if isempty(selection)
		S = stoc.S
		selection = [1:S;]
	end

    if haskey(options, :license)
		license = "distributed"
	else
		license = "new"
	end

	if haskey(options, :builtModel)
		allSubprobs = options[:builtModel]
		@assert isa(allSubprobs, Array{oneProblem})
		@assert length(allSubprobs) == length(selection)
		builtModel = true
		if isa(soln, solnType)
			pgstarter = soln.primal[:pg]
			hstarter = soln.primal[:h]
		elseif isa(soln, designType)
			pgstarter = soln.pg
			hstarter = soln.h
		else
			error("Unkown solution type.")
		end
	else
		builtModel = false
	end

	for s in selection
		if builtModel
			# This feature only works when all scnearios's feasibility are checked
			@assert length(selection) == stoc.S
			subprob = allSubprobs[s]
			for i = 1:(allSubprobs[s].param[:B])
				for j = 1:(allSubprobs[s].param[:T])
					setupperbound(allSubprobs[s].vars[:pg][i,j], pgstarter[i,j])
					setlowerbound(allSubprobs[s].vars[:pg][i,j], pgstarter[i,j])
					setupperbound(allSubprobs[s].vars[:h][i,j], hstarter[i,j])
					setlowerbound(allSubprobs[s].vars[:h][i,j], hstarter[i,j])
				end
			end
			# ============================= NEW MODIFIED SECTION ================================ #
			allSubprobs[s].param[:assDet] = ones(Int, param[:B], param[:T], 1)
			for i in 1:param[:B]
				for t in 1:param[:T]
					# Reconclude ass information
					expandVal = getupperbound(allSubprobs[s].vars[:pg][i,t])
					hardenVal = getupperbound(allSubprobs[s].vars[:h][i,t])

					if stoc.scenarios[s].data["SS"][i,t] - (param[:ProM][i]*hardenVal + param[:Ele][i]) > 0
						allSubprobs[s].param[:assDet][i,t,1] = 0
					end

					setupperbound(allSubprobs[s].vars[:a][i,t,1], param[:aslDet][i,t,s] * allSubprobs[s].param[:assDet][i,t,1])
					setlowerbound(allSubprobs[s].vars[:a][i,t,1], param[:aslDet][i,t,s] * allSubprobs[s].param[:assDet][i,t,1])
					setupperbound(allSubprobs[s].vars[:ap][i,t,1], param[:aslDet][i,t,s] * allSubprobs[s].param[:assDet][i,t,1] * expandVal)
					setlowerbound(allSubprobs[s].vars[:ap][i,t,1], param[:aslDet][i,t,s] * allSubprobs[s].param[:assDet][i,t,1] * expandVal)
					setupperbound(allSubprobs[s].vars[:ah][i,t,1], allSubprobs[s].param[:assDet][i,t,1] * hardenVal)
					setlowerbound(allSubprobs[s].vars[:ah][i,t,1], allSubprobs[s].param[:assDet][i,t,1] * hardenVal)
					setupperbound(allSubprobs[s].vars[:ass][i,t,1], allSubprobs[s].param[:assDet][i,t,1])
					setlowerbound(allSubprobs[s].vars[:ass][i,t,1], allSubprobs[s].param[:assDet][i,t,1])
				end
			end
			# ============================= NEW MODIFIED SECTION ================================ #
		else
			subprob = oneProblem()
			subprob = sbd_base_formulation(power, param, stoc, soln)
			subprob = attach_scenario(subprob, stoc, [s], exargs[:MODEL], 0.0, exargs, subprobType="tight")
			# ============================= NEW MODIFIED SECTION ================================ #
			subprob.param[:assDet] = ones(Int, param[:B], param[:T], 1)
			for i in 1:param[:B]
				for t in 1:param[:T]
					expandVal = getupperbound(subprob.vars[:pg][i,t])
					hardenVal = getupperbound(subprob.vars[:h][i,t])
					if stoc.scenarios[s].data["SS"][i,t] - (param[:ProM][i]*hardenVal + param[:Ele][i]) > 0
						subprob.param[:assDet][i,t,1] = 0
					end
					@assert	haskey(param, :aslDet)
					setupperbound(subprob.vars[:a][i,t,1], param[:aslDet][i,t,s] * subprob.param[:assDet][i,t,1])
					setlowerbound(subprob.vars[:a][i,t,1], param[:aslDet][i,t,s] * subprob.param[:assDet][i,t,1])
					setupperbound(subprob.vars[:ap][i,t,1], param[:aslDet][i,t,s] * subprob.param[:assDet][i,t,1] * expandVal)
					setlowerbound(subprob.vars[:ap][i,t,1], param[:aslDet][i,t,s] * subprob.param[:assDet][i,t,1] * expandVal)
					setupperbound(subprob.vars[:ah][i,t,1], subprob.param[:assDet][i,t,1] * hardenVal)
					setlowerbound(subprob.vars[:ah][i,t,1], subprob.param[:assDet][i,t,1] * hardenVal)
					setupperbound(subprob.vars[:ass][i,t,1], subprob.param[:assDet][i,t,1])
					setlowerbound(subprob.vars[:ass][i,t,1], subprob.param[:assDet][i,t,1])
				end
			end
			# ============================= NEW MODIFIED SECTION ================================ #
		end
		@objective(subprob.model, Min, 0)

		if config.PARALLEL
				config_solver(subprob.model, license=1, timelimit=20, showlog=0, focus="feasibility", presolve=1, threads=1)
		else
			config_solver(subprob.model, license=config.ENVS, timelimit=20, showlog=0, focus="feasibility", presolve=1, threads=1)
		end

		status = solve(subprob.model, suppress_warnings=true)

		if status == :Optimal
			feaCnt += 1
			push!(feaPool,s)
		else
			push!(infeaPool,s)
		end
	end

	ub_coverage = feaCnt / stoc.S
	# Reserve the output ub for later
	ub = get_design_cost(soln, param)

	return ub, ub_coverage, infeaPool, feaPool
end

"""
	Input ::
		power::Dict | param::Dict | stoc::stocType | soln | exargs::Dict | selection = [] | ...
	Usage ::
		Similar to check_feasible, this checks a certain design's slackness on all other scenarios. Get this function ready for sbd_norisk type of problem.
	Output ::
		UB | UB_Indicator | Infeasible Set | Feasible Set
"""
function check_slackness(power::Dict, param::Dict, stoc::stocType, soln, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)

	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"

	feaCnt = 0
	feaDict = Dict()
	infeaDict = Dict()
	maxSlack = -Inf
	maxSlackIdx = 0
	coverage = 0.0

	if isempty(selection) == true
		S = stoc.S
		selection = [1:S;]
	end

	slackPool = zeros(length(selection))

	if haskey(options, :builtModel)
		allSubprobs = options[:builtModel]
		@assert isa(allSubprobs, Array{oneProblem})
		@assert length(allSubprobs) >= length(selection) # All problems must be presented here
		builtModel = true
		if isa(soln, solnType)
			pgstarter = soln.primal[:pg]
			hstarter = soln.primal[:h]
		elseif isa(soln, designType)
			pgstarter = soln.pg
			hstarter = soln.h
		else
			error("Unkown solution type.")
		end
	else
		builtModel = false
	end

	# Performing feasibility check on all scenarios
	for s in 1:length(selection)
		if builtModel
			@assert length(selection) <= stoc.S
			slackProbScenario = pop!(selection)
			slackProb = allSubprobs[slackProbScenario]
			for i = 1:(allSubprobs[slackProbScenario].param[:B])
				for j = 1:(allSubprobs[slackProbScenario].param[:T])
					setupperbound(allSubprobs[slackProbScenario].vars[:pg][i,j], pgstarter[i,j])
					setlowerbound(allSubprobs[slackProbScenario].vars[:pg][i,j], pgstarter[i,j])
					setupperbound(allSubprobs[slackProbScenario].vars[:h][i,j], hstarter[i,j])
					setlowerbound(allSubprobs[slackProbScenario].vars[:h][i,j], hstarter[i,j])
				end
			end

			# ============================= NEW MODIFIED SECTION ================================ #
			allSubprobs[slackProbScenario].param[:assDet] = ones(Int, param[:B], param[:T], 1)
			for i in 1:param[:B]
				for t in 1:param[:T]
					# Reconclude ass information
					hardenVal = getupperbound(allSubprobs[slackProbScenario].vars[:h][i,t])
					expandVal = getupperbound(allSubprobs[slackProbScenario].vars[:pg][i,t])

					# Double checking the fixed first stage decision
					@assert expandVal == pgstarter[i,t]
					@assert hardenVal == hstarter[i,t]

					# Reconclude assDet parameter
					if stoc.scenarios[slackProbScenario].data["SS"][i,t] - (param[:ProM][i]*hardenVal + param[:Ele][i]) > 0
						allSubprobs[slackProbScenario].param[:assDet][i,t,1] = 0
					end

					@assert	haskey(param, :aslDet)
					setupperbound(allSubprobs[slackProbScenario].vars[:a][i,t,1], param[:aslDet][i,t,slackProbScenario]*allSubprobs[slackProbScenario].param[:assDet][i,t,1])
					setlowerbound(allSubprobs[slackProbScenario].vars[:a][i,t,1], param[:aslDet][i,t,slackProbScenario]*allSubprobs[slackProbScenario].param[:assDet][i,t,1])
					setupperbound(allSubprobs[slackProbScenario].vars[:ap][i,t,1], param[:aslDet][i,t,slackProbScenario]*allSubprobs[slackProbScenario].param[:assDet][i,t,1]*expandVal)
					setlowerbound(allSubprobs[slackProbScenario].vars[:ap][i,t,1], param[:aslDet][i,t,slackProbScenario]*allSubprobs[slackProbScenario].param[:assDet][i,t,1]*expandVal)
					setupperbound(allSubprobs[slackProbScenario].vars[:ah][i,t,1], allSubprobs[slackProbScenario].param[:assDet][i,t,1]*hardenVal)
					setlowerbound(allSubprobs[slackProbScenario].vars[:ah][i,t,1], allSubprobs[slackProbScenario].param[:assDet][i,t,1]*hardenVal)
					setupperbound(allSubprobs[slackProbScenario].vars[:ass][i,t,1], allSubprobs[slackProbScenario].param[:assDet][i,t,1])
					setlowerbound(allSubprobs[slackProbScenario].vars[:ass][i,t,1], allSubprobs[slackProbScenario].param[:assDet][i,t,1])
				end
			end
			# ============================= NEW MODIFIED SECTION ================================ #

		else
			slackProb = oneProblem()
			slackProb = sbd_base_formulation(power, param, stoc, soln)
			slackProbScenario = pop!(selection)
			slackProb = attach_scenario(slackProb, stoc, [slackProbScenario], exargs[:MODEL], 0.0, exargs, subprobType="slackness")

			# ============================= NEW MODIFIED SECTION ================================ #
			slackProb.param[:assDet] = ones(Int, param[:B], param[:T], 1)
			for i in 1:param[:B]
				for t in 1:param[:T]
					# Reconclude ass information
					expandVal = getupperbound(slackProb.vars[:pg][i,t])
					hardenVal = getupperbound(slackProb.vars[:h][i,t])

					# Reconclude assDet parameter
					if stoc.scenarios[slackProbScenario].data["SS"][i,t] - (param[:ProM][i]*hardenVal + param[:Ele][i]) > 0
						slackProb.param[:assDet][i,t,1] = 0
					end

					@assert	haskey(param, :aslDet)
					setupperbound(slackProb.vars[:a][i,t,1], param[:aslDet][i,t,slackProbScenario] * slackProb.param[:assDet][i,t,1])
					setlowerbound(slackProb.vars[:a][i,t,1], param[:aslDet][i,t,slackProbScenario] * slackProb.param[:assDet][i,t,1])
					setupperbound(slackProb.vars[:ap][i,t,1], param[:aslDet][i,t,slackProbScenario] * slackProb.param[:assDet][i,t,1] * expandVal)
					setlowerbound(slackProb.vars[:ap][i,t,1], param[:aslDet][i,t,slackProbScenario] * slackProb.param[:assDet][i,t,1] * expandVal)
					setupperbound(slackProb.vars[:ah][i,t,1], slackProb.param[:assDet][i,t,1] * hardenVal)
					setlowerbound(slackProb.vars[:ah][i,t,1], slackProb.param[:assDet][i,t,1] * hardenVal)
					setupperbound(slackProb.vars[:ass][i,t,1], slackProb.param[:assDet][i,t,1])
					setlowerbound(slackProb.vars[:ass][i,t,1], slackProb.param[:assDet][i,t,1])
				end
			end
			# ============================= NEW MODIFIED SECTION ================================ #
		end

		config_solver(slackProb.model, license=config.ENVS, timelimit=1800, showlog=0, focus="optimality", presolve=1, threads=config.WORKERTHREADS)

		status = solve(slackProb.model, suppress_warnings=true)
		if status == :Infeasible
			info("Original problem infeasible. Getting iis (this might take a while)...")
			error("ERROR|sbd.jl|check_slackness()|Formulation issue")
		end

		info("[>>>] Scenario $(slackProbScenario) with slackness $(getobjectivevalue(slackProb.model))")
		push!(slackPool, getobjectivevalue(slackProb.model))
		if abs(slackPool[end] - 0) > config.TOLERANCE
			infeaDict[slackProbScenario] = slackPool[end]
		else
			feaDict[slackProbScenario] = slackPool[end]
		end

		# Collect the "most infeaisble=" scenarios
		if slackPool[end] > maxSlack
			maxSlack = slackPool[end]
			maxSlackIdx = slackProbScenario
		end
	end

	return maxSlackIdx, maxSlack, infeaDict, feaDict
end

function sbd_master_formulation(power::Dict, param::Dict, stoc::stocType, exargs::Dict, overrideEps=exargs[:eps], cuts=[], prevSolution=[]; kwargs...)

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

"""
	Description goes here. This formulation doesn't fix solution.
"""
function sbd_subprob_formulation(prob::Dict, param::Dict, stoc::stocType, selection, exargs::Dict; kwargs...)

	options = Dict(kwargs)
	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"

	subprob = oneProblem()
	subprob = sbd_base_formulation(prob, param, stoc)
	subprob = attach_scenario(subprob, stoc, selection, exargs[:MODEL], 0.0, exargs, subprobType=subprobType)
	subprob = trim_bounds(subprob, stoc, exargs, selection=selection)

    return subprob
end
