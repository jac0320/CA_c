################################################################################
#	This is a specific problem :: Adaptive Incremental Design                  #
################################################################################

"""
	This is the basic formulation mainly utilized for solve the problem deterministically.
	Inputs -> power::Dict | param::Dict | stoc::stocType | exargs::Dict | ...
"""
function base_formulation(prob::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

	base = oneProblem()
	m = init_model_solver()
	param = check_parameter_intactness(param)

	# Attach parameters to the oneProblem structure
	base.param = param
	base.samples = stoc

	# Flip out some commonly used dimensions
	S = base.param[:S] = stoc.S # Be careful with this
	T = base.param[:T] = stoc.T
	B = base.param[:B]
	L = base.param[:L]

	MDouble = config.BIGMDOUBLE;
	MInt = config.BIGMINT

	# ========================== Record Model Info ========================== #
	base.name = prob["name"]
	base.stage = 1
	base.T = stoc.T
	base.S = stoc.S
	base.vars = Dict()
	# =========================== Setup the Model =========================== #

	# Variables :: ========================================================== #
	# Introduce the necessary variables
	base.vars[:pg] = @variable(m, 0<=pg[i=1:B, 1:T]<=param[:Pgbar][i], Int)
	base.vars[:h] = @variable(m, 0<=h[i=1:B, 1:T]<=param[:Hbar][i], Int)
	base.vars[:f] = @variable(m, f[1:T, 1:S], Bin)
	base.vars[:a] = @variable(m, a[1:B, 1:T, 1:S], Bin)

	base.vars[:ass] = @variable(m, ass[1:B, 1:T, 1:S], Bin)
	base.vars[:ap] = @variable(m, ap[1:B, 1:T, 1:S], Int)
	base.vars[:ah] = @variable(m, ah[1:B, 1:T, 1:S], Int)

	# Constraints ============================================================ #
	# (First Stage) Incremental Design Geneartion || Hardning
	@constraint(m, [i=1:B], pg[i,1] >= param[:Pg0][i]);
	@constraint(m, [i=1:B], h[i,1] >= param[:H0][i]);
	@constraint(m, [i=1:B,t=2:T], pg[i,t-1] <= pg[i,t]);
	@constraint(m, [i=1:B,t=2:T], h[i,t-1] <= h[i,t]);

	# Above are reduced due to the information can be determinied deterministically
	@constraint(m, [i=1:B, t=1:T, s=1:S],
		stoc.scenarios[s].data["SS"][i,t] * (2*base.vars[:ass][i,t,s]-1) - 2*base.param[:Ele][i]*base.vars[:ass][i,t,s] - 2*base.param[:ProM][i]*base.vars[:ah][i,t,s] + base.param[:Ele][i] + base.param[:ProM][i]*base.vars[:h][i,t] <= 0)

	# McCormick Relaxation => ap = a * pg # Replace this moduleswith methods in relax.jl
	for i in 1:B
		for t in 1:T
			for s in 1:S
				@constraint(m, a[i,t,s] == base.param[:aslDet][i,t,s] * base.vars[:ass][i,t,s])
				mccormick(m, base.vars[:ap][i,t,s], base.vars[:a][i,t,s], base.vars[:pg][i,t], 0, 1, param[:Pg0][i], param[:Pgbar][i]);
				mccormick(m, base.vars[:ah][i,t,s] ,base.vars[:ass][i,t,s], base.vars[:h][i,t], 0, 1, param[:H0][i], param[:Hbar][i]);
			end
		end
	end

	@objective(m, Min, sum(param[:Cg][i,1]*(pg[i,1]-param[:Pg0][i]) for i=1:B)
						+ sum(param[:Ch][i,1]*(h[i,1]-param[:H0][i]) for i=1:B)
	                    + sum(param[:Cg][i,t]*(pg[i,t]-pg[i,t-1]) for i=1:B, t=2:T)
						+ sum(param[:Ch][i,t]*(h[i,t]-h[i,t-1]) for i=1:B, t=2:T));

	# ======================================================================== #
	base.model = m;
	return base
end

"""
	Upon the base_formulation(), this functions adds capacity characteristic to the problem.
	Inputs -> power::Dict | exargs::Dict | ...
"""
function capacity_characteristic(prob::oneProblem, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)

	if haskey(options, :subprobType)
		subprobType = options[:subprobType]
	else
		subprobType = "free"
	end

	B = prob.param[:B]
	T = prob.T
	S = prob.S

	if isempty(selection)
		selection = [1:S;]
	end

	@assert haskey(prob.vars, :ap)

	MInt = config.BIGMINT
	MDouble = config.BIGMDOUBLE

	prob.vars[:pdv] = @variable(prob.model, pdv[1:B, 1:T, 1:S]>=0)

	if subprobType == "free"
		@assert haskey(prob.vars, :f)
		@assert haskey(prob.vars, :pdv)

		@constraint(prob.model, [t=1:T,s=1:S],
			sum(prob.param[:PgUB][i] * prob.vars[:ap][i,t,s] for i=1:B) - sum(prob.vars[:pdv][i,t,s] for i=1:B) >= -sum(prob.param[:aslDetPd][i,t,s] for i=1:B)*(1-prob.vars[:f][t,s]))
		@constraint(prob.model, [t=1:T,s=1:S],
			sum(prob.param[:aslDetPd][i,t,s] for i=1:B) * prob.vars[:f][t,s] >= sum(prob.param[:PgUB][i]*prob.vars[:ap][i,t,s] for i=1:B) - sum(prob.vars[:pdv][i,t,s] for i=1:B))

	elseif subprobType == "tight"
		@constraint(prob.model, [t=1:T,s=1:S],
			sum(prob.param[:PgUB][i] * prob.vars[:ap][i,t,s] for i=1:B) >= sum(prob.vars[:pdv][i,t,s] for i=1:B))

	elseif subprobType == "slackness"
		prob.vars[:scapslack] = @variable(prob.model, capslack[1:T,1:S]>=0)
		@constraint(prob.model, [t=1:T,s=1:S],
			sum(prob.param[:PgUB][i] * prob.vars[:ap][i,t,s] for i=1:B) + prob.vars[:capslack][t,s] >= sum(prob.param[:Pd][i,t] for i=1:B))

		@objective(prob.model, Min, sum(prob.vars[:capslack]))

	else
		error("ERROR|general.jl|capacity_characteristic()|Unknown subproblem type.")
	end

	return prob
end

"""
	Upon the base_formulation(), this function adds network characteristic to the probelm.
	Such characteristic are with binary indicators for feasibility indication.
	Inputs -> power::Dict | exargs::Dict | ...
"""
function network_characteristic(prob::oneProblem, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)

	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	haskey(options, :logical) ? logical = options[:logical] : logical = true
	(isempty(selection)) && (selection = [1:prob.param[:S];])

	# Commonly used parameters :: keep an eye on their size
	B = prob.param[:B]
	T = prob.T
	S = length(selection)
	L = prob.param[:L]

	@assert prob.S == prob.param[:S] # Delete this later
	@assert prob.S == length(selection)
	@assert haskey(prob.vars, :ap)
	@assert haskey(prob.vars, :ass)

	MInt = config.BIGMINT
	MDouble = config.BIGMDOUBLE

	#Power Flow i->j | Actual Generation | Flow Balance Feasible Indicator

	# If a bus is flooded, then it linked thermal limits will be cutted off
	# If a bus is undersea permentaly, also cut off the associated thermal limits
	prob.vars[:p] = @variable(prob.model, p[i=1:B,j=1:B,1:T,1:S; prob.param[:EDGE][i,j]==1])

	for i in 1:B
		# Line switch off given sea level rise
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["in"]],
			prob.vars[:p][j,i,t,s] <= prob.vars[:ass][i,t,s] * prob.param[:Lcap][j,i])
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["in"]],
			prob.vars[:p][j,i,t,s] >= prob.vars[:ass][i,t,s] * (-prob.param[:Lcap][j,i]))
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["out"]],
			prob.vars[:p][i,j,t,s] <= prob.vars[:ass][i,t,s] * prob.param[:Lcap][i,j])
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["out"]],
			prob.vars[:p][i,j,t,s] >= prob.vars[:ass][i,t,s] * (-prob.param[:Lcap][i,j]))
		# Line swtich off given storm surge
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["in"]],
			prob.vars[:p][j,i,t,s] <= prob.param[:aslDet][i,t,selection[s]] * prob.param[:Lcap][j,i])
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["in"]],
			prob.vars[:p][j,i,t,s] >= prob.param[:aslDet][i,t,selection[s]] * (-prob.param[:Lcap][j,i]))
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["out"]],
			prob.vars[:p][i,j,t,s] <= prob.param[:aslDet][i,t,selection[s]] * prob.param[:Lcap][i,j])
		@constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in prob.param[:Edge][i]["out"]],
			prob.vars[:p][i,j,t,s] >= prob.param[:aslDet][i,t,selection[s]] * (-prob.param[:Lcap][i,j]))
	end

	# Generation applicability indicator logics
	prob.vars[:pdv] = @variable(prob.model, pdv[1:B, 1:T, 1:S]>=0) # Variable Demand
	prob.vars[:psg] = @variable(prob.model, psg[1:B, 1:T, 1:S]>=0) # Variable generation output

	# Varibale generation ouput upper bounded by total capaciated units
	@constraint(prob.model, [i=1:B,t=1:T,s=1:S],
		prob.vars[:psg][i,t,s] <= prob.param[:PgUB][i] * prob.vars[:ap][i,t,s])

	if subprobType == "free"

		@assert haskey(prob.vars, :f)

		prob.vars[:fn] = @variable(prob.model, fn[1:B, 1:T, 1:S], Bin)

		# Flow Balance Constraints
		for i in 1:B
			@constraint(prob.model, [t=1:T,s=1:S],
				sum(prob.vars[:p][j,i,t,s] for j in prob.param[:Edge][i]["in"])- sum(prob.vars[:p][i,j,t,s] for j in prob.param[:Edge][i]["out"]) + prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s] <= MInt * prob.vars[:fn][i,t,s])
			@constraint(prob.model, [t=1:T,s=1:S],
				-(sum(prob.vars[:p][j,i,t,s] for j in prob.param[:Edge][i]["in"]) - sum(prob.vars[:p][i,j,t,s] for j in prob.param[:Edge][i]["out"]) + prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s]) <= MInt * prob.vars[:fn][i,t,s])
		end

		if logical	# If this is skipped, then it is a DC network, similar logical constraints will be added in dc part
			@constraint(prob.model, [s=1:S,t=1:T],
				B*(prob.vars[:f][t,s]-1) <= sum(prob.vars[:fn][i,t,s] for i=1:B))
			@constraint(prob.model, [s=1:S,t=1:T],
				sum(prob.vars[:fn][i,t,s] for i=1:B) <= B*(1-prob.vars[:f][t,s]))
		end

	elseif subprobType == "tight"

		for i in 1:B
			@constraint(prob.model, [t=1:T, s=1:S],
				sum(prob.vars[:p][j,i,t,s] for j in prob.param[:Edge][i]["in"]) - sum(prob.vars[:p][i,j,t,s] for j in prob.param[:Edge][i]["out"]) + prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s] == 0)
		end

	elseif subprobType == "slackness"

		prob.vars[:flowSlackpos] = @variable(prob.model, flowSlackpos[1:B,1:T,1:S]>=0)
		prob.vars[:flowSlackneg] = @variable(prob.model, flowSlackneg[1:B,1:T,1:S]>=0)

		for i in 1:B
			@constraint(prob.model, [t=1:T,s=1:S],
				sum(prob.vars[:p][j,i,t,s] for j in prob.param[:Edge][i]["in"]) - sum(prob.vars[:p][i,j,t,s] for j in prob.param[:Edge][i]["out"]) + prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s] + prob.vars[:flowSlackpos][i,t,s] - prob.vars[:flowSlackneg][i,t,s] == 0)
		end

		@objective(prob.model, Min,
			 sum(prob.vars[:flowSlackpos]) + sum(prob.vars[:flowSlackneg]))

	else
		error("ERROR|general.jl|network_characteristic()|Unkown subproblem type.")
	end

	return prob
end

"""
	Upon the base_formulation(), this function adds network characteristic to the probelm.
	Such characteristic are with binary indicators for feasibility indication.
	Inputs -> power::Dict | exargs::Dict | ...

	Higher level
"""
function dc_characteristic(prob::oneProblem, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)

	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"

	# Extra dimension from parameter dicstionary
	B = prob.param[:B]
	L = prob.param[:L]

	# Consistence checking
	@assert prob.param[:T] == prob.T
	@assert prob.param[:S] == prob.S

	T = prob.T
	S = prob.S
	E = length(prob.param[:line])

	(isempty(selection)) && (selection = [1:S;])

	# populate lineX value
	lineX = zeros(Float64, B, B)
	for l in 1:E
		lineX[prob.param[:line][l]["f_bus"],prob.param[:line][l]["t_bus"]] = prob.param[:line][l]["br_x"]
	end
	prob.param[:lineX] = lineX / 100

	@assert haskey(prob.param, :RefBus)

	prob = network_characteristic(prob, exargs, selection, subprobType=subprobType, logical=false)

	# Fetch symbols from the established model
	@assert haskey(prob.vars, :p)
	@assert haskey(prob.vars, :f)

	# Add Angles Variable
	prob.vars[:theta] = @variable(prob.model,
		prob.param[:AngleLimit]>=theta[1:B, 1:T, 1:S]>=-prob.param[:AngleLimit])

	# Reference bus with zero angle
	@constraint(prob.model, [s=1:S, t=1:T],
		prob.vars[:theta][prob.param[:RefBus],t,s] == 0)

	# Voltage are assumed to be one; assumptions, parameters never change during all time periods
	# Flow Balance Equation added p[(i,j),s,t] = (theta[i]-theta[j]) / x[(i,j)]
	# There is room for improvement
	if subprobType == "free"

		# Make sure another logical constraint is in the model
		@assert haskey(prob.vars, :fn)

		prob.vars[:fln] = @variable(prob.model, fln[i=1:B,j=1:B,1:T,1:S; prob.param[:EDGE][i,j] == 1], Bin)
		prob.vars[:assss] = @variable(prob.model, assss[i=1:B,j=1:B,1:T,1:S; prob.param[:EDGE][i,j] == 1], Bin)

		@constraint(prob.model, [t=1:T,s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1], # Seal Upper bound
			prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/prob.param[:lineX][i,j] <=
			(prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * prob.vars[:fln][i,j,t,s] + 100*(prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] * prob.param[:aslDet][i,t,selection[s]] * prob.param[:aslDet][j,t,selection[s]]))

		@constraint(prob.model, [t=1:T,s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1], # Seal lower bound
			-(prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/prob.param[:lineX][i,j]) <=
			  (prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * prob.vars[:fln][i,j,t,s] + 100*(prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] * prob.param[:aslDet][i,t,selection[s]] * prob.param[:aslDet][j,t,selection[s]]))

	  	@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1],
			prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] <= prob.param[:AngleShiftLimit])
		@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1],
			prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] >= -prob.param[:AngleShiftLimit])


		# A huge extra system of constraints if DC model is included
		for i in 1:B
			for j in 1:B
				if prob.param[:EDGE][i,j] == 1
					if "surge-load-shed" in exargs[:FEATURES]
						for t in 1:T
							for s in 1:S
								mccormick(prob.model, prob.vars[:assss][i,j,t,s], prob.vars[:ass][i,t,s], prob.vars[:ass][j,t,s], 0, 1, 0, 1)
							end
						end
					else
						error("Must have surge-load-shed formulation.")
					end
				end
			end
		end

		# Add Logical Feasibility: this means the logical constraints in network functin is passed
		@constraint(prob.model, [s=1:S,t=1:T],
			(B+L)*(prob.vars[:f][t,s]-1) <= sum(prob.vars[:fn][i,t,s] for i=1:B) + sum(prob.vars[:fln][i,j,t,s] for i=1:B,j=1:B if prob.param[:EDGE][i,j]==1))
		@constraint(prob.model, [s=1:S,t=1:T],
			sum(prob.vars[:fn][i,t,s] for i=1:B) + sum(prob.vars[:fln][i,j,t,s] for i=1:B, j=1:B if prob.param[:EDGE][i,j]==1) <= (B+L)*(1-prob.vars[:f][t,s]))

	elseif subprobType == "tight"

		prob.vars[:assss] = @variable(prob.model, assss[i=1:B,j=1:B,1:T,1:S; prob.param[:EDGE][i,j] == 1], Bin)

		# ======================= Old Implementation ======================= #
		@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1], # Seal Upper bound
			prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/prob.param[:lineX][i,j] <= 100*(prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] * prob.param[:aslDet][i,t,selection[s]] * prob.param[:aslDet][j,t,selection[s]]) )

		@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1], # Seal lower bound
		   -(prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/prob.param[:lineX][i,j]) <= 100*(prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] * prob.param[:aslDet][i,t,selection[s]] * prob.param[:aslDet][j,t,selection[s]]))
		# ======================= Old Implementation ======================= #

		@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1],
			prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] <= prob.param[:AngleShiftLimit])
		@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1],
			prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] >= -prob.param[:AngleShiftLimit])

		# A huge extra system of constraints if DC model is included
		for i in 1:B
			for j in 1:B
				if prob.param[:EDGE][i,j] == 1
					for t in 1:T
						for s in 1:S
							mccormick(prob.model, prob.vars[:assss][i,j,t,s], prob.vars[:ass][i,t,s], prob.vars[:ass][j,t,s], 0, 1, 0, 1)
						end
					end
				end
			end
		end


	elseif subprobType == "slackness"

		prob.vars[:flowEquSlackpos] = @variable(prob.model,
			flowEquSlackpos[i=1:B,j=1:B,t=1:T,s=1:S; prob.param[:EDGE][i,j] == 1]>=0)
		prob.vars[:flowEquSlackneg] = @variable(prob.model,
			flowEquSlackneg[i=1:B,j=1:B,t=1:T,s=1:S; prob.param[:EDGE][i,j] == 1]>=0)

		 # Setup the production of both ass
		prob.vars[:assss] = @variable(prob.model, assss[i=1:B,j=1:B,1:T,1:S; prob.param[:EDGE][i,j]==1], Bin)

		@constraint(prob.model, [t=1:T,s=1:S,i=1:B,j=1:B; prob.param[:EDGE][i,j]==1], # Seal Upper bound
		 	prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/prob.param[:lineX][i,j] + prob.vars[:flowEquSlackpos][i,j,t,s] - prob.vars[:flowEquSlackneg][i,j,t,s] <=
			(prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] * prob.param[:aslDet][i,t,selection[s]] * prob.param[:aslDet][j,t,selection[s]]) )

		@constraint(prob.model, [t=1:T,s=1:S,i=1:B,j=1:B; prob.param[:EDGE][i,j]==1], # Seal lower bound
		 	-(prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/prob.param[:lineX][i,j] + prob.vars[:flowEquSlackpos][i,j,t,s] - prob.vars[:flowEquSlackneg][i,j,t,s]) <=
			(prob.param[:Lcap][i,j] + prob.param[:AngleShiftLimit] / prob.param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] * prob.param[:aslDet][i,t,selection[s]] * prob.param[:aslDet][j,t,selection[s]]) )

	  	@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1],
  			prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] <= prob.param[:AngleShiftLimit])
  		@constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; prob.param[:EDGE][i,j] == 1],
  			prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] >= -prob.param[:AngleShiftLimit])

		@objective(prob.model, Min, sum(prob.vars[:flowSlackneg]) + sum(prob.vars[:flowSlackpos]) + sum(prob.vars[:flowEquSlackneg]) + sum(prob.vars[:flowEquSlackpos]))

		# A huge extra system of constraints if DC model is included
		for i in 1:B
			for j in 1:B
				if prob.param[:EDGE][i,j] == 1
					for t in 1:T
						for s in 1:S
							mccormick(prob.model, prob.vars[:assss][i,j,t,s], prob.vars[:ass][i,t,s], prob.vars[:ass][j,t,s], 0, 1, 0, 1)
						end
					end
				end
			end
		end

	else
		error("ERROR|general.jl|dc_characteristic()|Unkown subprob type.")
	end

	return prob
end

"""
	This function generate extra constraints to the problem to address some available features upon the established
	formulation
"""
function formulation_features(prob::oneProblem, stoc::stocType, exargs::Dict, overrideEps=exargs[:eps], selection=[]; kwargs...)

	options = Dict(kwargs)

	S = prob.param[:S]
	B = prob.param[:B]
	T = prob.param[:T]
	eps = overrideEps

	if isempty(selection)
		selection = [1:S;]
	end

	if haskey(options, :subprobType)
		subprobType = options[:subprobType]
	else
		subprobType = "free"
	end

	MInt = config.BIGMINT
	MDouble = config.BIGMDOUBLE

	# Supplyment choice
	feaLoadShift = false

	for feature in exargs[:FEATURES]
		if feature == "sample-based-risk"
			@assert haskey(prob.vars, :f)
			prob.vars[:fs] = @variable(prob.model,fs[1:S],Bin)
			@constraint(prob.model, [s=1:S],
				T * prob.vars[:fs][s] <= sum(prob.vars[:f][t,s] for t=1:T))
			@constraint(prob.model, risk,
				sum(fs[s] for s=1:S) >= ceil(S * (1 - eps)))

		elseif feature == "time-based-risk"
			@assert haskey(prob.vars, :f)
			@constraint(prob.model, [t=1:T],
				sum(prob.vars[:f][t,s] for s=1:S) >= S*(1-eps))

		elseif feature == "short-term-budget"
			@assert haskey(prob.vars, :pg)
			@assert haskey(prob.vars, :h)
			@assert haskey(prob.param, :budget)

			println("Budget parameter used");
			@constraint(prob.model,
				sum(prob.param[:Cg][i,1]*(prob.vars[:pg][i,1]-prob.param[:Pg0][i]) for i=1:B) + sum(prob.param[:Ch][i,1]*(prob.vars[:h][i,1]-prob.param[:H0][i]) for i=1:B) <= prob.param[:budget])
			@constraint(prob.model, [t=2:T],
				sum(prob.param[:Cg][i]*(prob.vars[:pg][i,t]-prob.vars[:pg][i,t-1]) for i=1:B) + sum(prob.param[:Ch][i,t]*(prob.vars[:h][i,t]-prob.vars[:h][i,t-1]) for i=1:B) <= prob.param[:budeget])

		elseif feature == "long-term-budget"
			@assert haskey(prob.vars, :pg)
			@assert haskey(prob.vars, :h)
			@assert haskey(prob.param, :budget)
			@constraint(prob.model,
				sum(param[:Cg][i,1]*(prob.vars[:pg][i,1]-param[:Pg0][i]) for i=1:B) + sum(param[:Ch][i,1]*(prob.vars[:h][i,1]-param[:H0][i]) for i=1:B) + sum(param[:Cg][i]*(prob.vars[:pg][i,t]-prob.vars[:pg][i,t-1]) for i=1:B,t=2:T) + sum(param[:Ch][i,t]*(prob.vars[:h][i,t]-prob.vars[:h][i,t-1]) for i=1:B,t=2:T) <= prob.param[:budget]);

		elseif feature == "common-hit-policy"

			bigMInt = config.BIGMINT

			# Formulation follows section §7.3
			@assert haskey(prob.vars, :ass)
			@assert haskey(prob.vars, :pg)

			prob.vars[:pgb] = @variable(prob.model, pgb[1:B, 1:T], Bin)
			prob.vars[:v] = @variable(prob.model, v>=0)

			@constraint(prob.model, [i=1:B, t=1:T],
				sum((1-prob.vars[:ass][i,t,s]) for s=1:S)-v <= bigMInt*(1-prob.vars[:pgb][i,t]))
			@constraint(prob.model, [i=1:B, t=1:T],
				sum((1-prob.vars[:ass][i,t,s]) for s=1:S)-v >= -bigMInt*prob.vars[:pgb][i,t])
			@constraint(prob.model, [i=1:B],
				prob.vars[:pg][i,1] <= prob.param[:Pg0] + prob.param[:PgBar]*prob.vars[:pgb][i,1])
			@constraint(prob.model, [i=1:B, t=2:T],
				prob.vars[:pg][i,t] <= prob.vars[:pg][i,t-1] + prob.param[:PgBar]*prob.vars[:pgb][i,t])

		elseif feature == "surge-load-shed"

			if exargs[:MODEL] != capacity_characteristic

				feaLoadShift = true
				if subprobType == "free"
					@assert haskey(prob.vars, :f)
					@constraint(prob.model, [t=1:T, s=1:S],
				 		sum(prob.vars[:pdv][i,t,s] for i=1:B) >= prob.param[:SHEDLambda] * sum(prob.param[:aslDetPd][i,t,selection[s]] for i=1:B) * prob.vars[:fs][s])
					@constraint(prob.model, dispLb[i=1:B, t=1:T, s=1:S],
						prob.vars[:pdv][i,t,s] <= prob.param[:aslDetPd][i,t,selection[s]] * prob.vars[:ass][i,t,s])

				elseif subprobType == "slackness"
					prob.vars[:slsSlack] = @variable(prob.model, slsslack[1:T,1:S]>=0)
					@constraint(prob.model, shed[t=1:T, s=1:S],
						sum(prob.vars[:pdv][i,t,s] for i=1:B) + prob.vars[:slsSlack][t,s] >= prob.param[:SHEDLambda] * sum(prob.param[:aslDetPd][i,t,selection[s]] for i=1:B))
					@constraint(prob.model, dispLb[i=1:B, t=1:T, s=1:S],
						prob.vars[:pdv][i,t,s] <= prob.param[:aslDetPd][i,t,selection[s]] * prob.vars[:ass][i,t,s])

					if exargs[:MODEL] == dc_characteristic
						@objective(prob.model, Min,
							sum(prob.vars[:slsSlack]) + sum(prob.vars[:flowSlackpos]) + sum(prob.vars[:flowSlackneg]) +  sum(prob.vars[:flowEquSlackpos]) + sum(prob.vars[:flowEquSlackneg]))
					elseif exargs[:MODEL] == network_characteristic
						@objective(prob.model, Min,
							sum(prob.vars[:slsSlack]) + sum(prob.vars[:flowSlackpos]) + sum(prob.vars[:flowSlackneg]))
					else
						error("Unkown model type for resettng objective in surge-load-shed condition.")
					end

				else subprobType == "tight"
					@constraint(prob.model, [i=1:B, t=1:T, s=1:S],
						prob.vars[:pdv][i,t,s] <= prob.param[:aslDetPd][i,t,selection[s]] * prob.vars[:ass][i,t,s])
					@constraint(prob.model, [t=1:T, s=1:S],
						sum(prob.vars[:pdv][i,t,s] for i=1:B) >= prob.param[:SHEDLambda] * sum(prob.param[:aslDetPd][i,t,selection[s]] for i=1:B))
				end

			else
				@constraint(prob.model, [t=1:T, s=1:S],
					sum(prob.vars[:pdv][i,t,s] for i=1:B) >= prob.param[:SHEDLambda]*sum(prob.param[:aslDetPd][i,t,s] for i=1:B))
				@constraint(prob.model, dispLb[i=1:B, t=1:T, s=1:S],
					prob.vars[:pdv][i,t,s] <= prob.param[:aslDetPd][i,t,selection[s]]*prob.vars[:ass][i,t,s])
			end
		end

	end

	return prob
end
