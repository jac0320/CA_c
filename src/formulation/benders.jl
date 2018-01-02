"""
    Builds the add variables to the model and returns variable set.
"""
function climate_benders_master_variables(m::JuMP.Model; kwargs...)

	data = Dict(kwargs)
	param = data[:package][:param]

    B = param[:B]
    T = param[:T]
    S = param[:S]

    Xvars = Dict()

    Xvars[:obj] = @variable(m, obj>=0)	# First Stage Cost
    Xvars[:pg] 	= @variable(m, 0<=pg[i=1:B, 1:T]<=param[:Pgbar][i], Int)
    Xvars[:h] 	= @variable(m, 0<=h[i=1:B, 1:T]<=param[:Hbar][i], Int)
    Xvars[:f]	= @variable(m, f[1:S], Bin)
    Xvars[:a]	= @variable(m, a[1:B, 1:T, 1:S], Bin)
    # This variable is no longer needed, but kept here
    # Xvars[:asl]	= @variable(m, asl[i=1:B, t=1:T, s=1:S] == param[:aslDet][i,t,s])
    Xvars[:ass]	= @variable(m, ass[1:B, 1:T, 1:S], Bin)
    Xvars[:ah]  = @variable(m, ah[1:B, 1:T, 1:S], Int)
    Xvars[:eta] = @variable(m, eta>=0)	#Lower Bound Coded	# Expected Recourse Function

    return m, Xvars
end

function climate_benders_master_constraints(m::JuMP.Model, Xvars::Dict; kwargs...)

    data = Dict(kwargs)

    # Get the frequently used info
    param = data[:package][:param]
	eps = data[:package][:exargs][:eps]

    B = param[:B]
    T = param[:T]
    S = data[:package][:stoc].S

    # (First Stage) Incremental Design Geneartion || Hardning
    @constraint(m, [i=1:B],
        Xvars[:pg][i,1] >= param[:Pg0][i])
    @constraint(m, [i=1:B],
        Xvars[:h][i,1] >= param[:H0][i])
    @constraint(m, [i=1:B, t=2:T],
        Xvars[:pg][i,t-1] <= Xvars[:pg][i,t])
    @constraint(m, [i=1:B, t=2:T],
        Xvars[:h][i,t-1] <= Xvars[:h][i,t])

    # (First Stage) Risk constraint (Time-based-risk)
    @constraint(m,
        sum(Xvars[:f][s] for s=1:S) >= S * (1-eps))

    # Feasibility Indicator Rules
    # @constraint(m, [i=1:B, t=1:T, s=1:S],
    #     Xvars[:a][i,t,s] >= Xvars[:asl][i,t,s] + Xvars[:ass][i,t,s] - 1)
    # @constraint(m, [i=1:B, t=1:T, s=1:S],
    #     Xvars[:a][i,t,s] <= Xvars[:asl][i,t,s])
    # @constraint(m, [i=1:B, t=1:T, s=1:S],
    #     Xvars[:a][i,t,s] <= Xvars[:ass][i,t,s])
    #  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #    Replaced by given asl is parameters now
    #  ___________________________________________
	# Getting rid of the variable asl => it can be deterministically determined
	@constraint(m, [i=1:B, t=1:T, s=1:S],
		Xvars[:a][i,t,s] == param[:aslDet][i,t,s] * Xvars[:ass][i,t,s])

    # TODO :: make this cleaner
    for i in 1:B
		for t in 1:T
			for s in 1:S
				mccormick(m, Xvars[:ah][i,t,s] ,Xvars[:ass][i,t,s], Xvars[:h][i,t], 0, 1, param[:H0][i], param[:Hbar][i])
			end
		end
	end

    # Design Cost Equation
    @constraint(m,
        Xvars[:obj] == sum(param[:Cg][i,1]*(Xvars[:pg][i,1]-param[:Pg0][i]) for i=1:B) + sum(param[:Ch][i,1]*(Xvars[:h][i,1]-param[:H0][i]) for i=1:B) + sum(param[:Cg][i,t]*(Xvars[:pg][i,t]-Xvars[:pg][i,t-1]) for i=1:B,t=2:T) + sum(param[:Ch][i,t]*(Xvars[:h][i,t]-Xvars[:h][i,t-1]) for i=1:B,t=2:T))

    return m
end

function climate_bender_master_objective(m::JuMP.Model, Xvars::Dict; kwargs...)

    options = Dict(kwargs)

    if options[:hasEta]
        @objective(m, Min, Xvars[:obj] + Xvars[:eta])
    else
        @objective(m, Min, Xvars[:obj])
    end

    return m

end

function climate_benders_subprob_variables(m::JuMP.Model; kwargs...)

	data = Dict(kwargs)

	param = data[:package][:param]

	B = param[:B]
	T = param[:T]

	Yvars = Dict()

	# Introduce the necessary variables
	# Power Flow i->j | Actual Generation | Flow Balance Feasible Indicator
	Yvars[:p] 	= @variable(m, -param[:Lcap][i,j]<=p[i=1:B,j=1:B,1:T; param[:EDGE][i,j]==1]<=param[:Lcap][i,j]);
	Yvars[:psg] = @variable(m, psg[1:B, 1:T]>=0)
	Yvars[:pdv] = @variable(m, pdv[1:B, 1:T]>=0)
	Yvars[:ap]	= @variable(m, ap[1:B, 1:T]>=0)		# generation applicable indicator
	Yvars[:etap] = @variable(m, etap[1:B, 1:T]>=0)
	Yvars[:etan] = @variable(m, etan[1:B, 1:T]>=0)

	return m, Yvars

end

function climate_benders_subprob_constraints(m::JuMP.Model, Yvars::Dict, Xvars::Dict, scenario::scenarioType; kwargs...)

	data = Dict(kwargs)

    # Get the frequently used info :: hopefully julia won't copy it
    param = data[:package][:param]
	eps = data[:package][:exargs][:eps]

    B = param[:B]
    T = param[:T]
	M = config.BIGMDOUBLE
	s = scenario.ind

	# These reference can be smarter when one randomness showed up in one constraint
	stoMap = Dict()
	stoMap[:T] = []
	stoMap[:h] = []
	stoMap[:q] = []

	# =============== Constraints =============== #
	# Operation indicators
    @constraint(m, [i=1:B, t=1:T],
        Yvars[:ap][i,t] <= Xvars[:a][i,t,s] * param[:Pgbar][i])
    @constraint(m, [i=1:B, t=1:T],
        Yvars[:ap][i,t] <= Xvars[:pg][i,t])
    @constraint(m, [i=1:B ,t=1:T],
        -Yvars[:ap][i,t] <= -Xvars[:pg][i,t] - (Xvars[:a][i,t,s] - 1) * param[:Pgbar][i])

	@constraint(m, ssRef[i=1:B, t=1:T],
        scenario.data["SS"][i,t]*(2*Xvars[:ass][i,t,s]-1) - 2*param[:Ele][i]*Xvars[:ass][i,t,s] - 2*param[:ProM][i]*Xvars[:ah][i,t,s] + param[:Ele][i] + param[:ProM][i]*Xvars[:h][i,t] <= 0)

	# TODO : This is a very ugly set of code, but it solve the problem
	for i=1:B
		for t=1:T
			# Each item records [random symbol, associated element, coefficients]
			# Feed these into a data structure
            # Be very careful !!!!
            # notice that T transfer matrix should be kept with original sige since the negative is taken later
			push!(stoMap[:T], [["SS",i,t], ssRef[i,t], Xvars[:ass][i,t,s], 2])
			push!(stoMap[:h], [["SS",i,t], ssRef[i,t].idx, 1])
		end
	end

    @constraint(m, [i=1:B, t=1:T],
        sum(Yvars[:p][j,i,t] for j in param[:Edge][i]["in"]) - sum(Yvars[:p][i,j,t] for j in param[:Edge][i]["out"]) + Yvars[:psg][i,t] - Yvars[:pdv][i,t] - Yvars[:etan][i,t] + Yvars[:etap][i,t] == 0)

	@constraint(m, [i=1:B, t=1:T],
        Yvars[:etan][i,t] <= M * (1 - Xvars[:f][s]))

    @constraint(m, [i=1:B, t=1:T],
        Yvars[:etap][i,t] <= M * (1 - Xvars[:f][s]))

	@constraint(m, [i=1:B, t=1:T],
        Yvars[:psg][i,t] <= param[:PgUB][i]*Yvars[:ap][i,t])

	# Determining demand pdv TODO: construct the same thing here for demand constraints
	for i in 1:B
    	SHIFT = zeros(Float64, T)
        for t in 1:T
			nom = 0.0
			denom = 1
			for j in 1:B
				# All loss load due to SLR
				nom += param[:Pd][j,t] * (1-param[:aslDet][j, t, s])
				# Current implementation indicates a uniform shift
				denom += param[:aslDet][j,t,s]
			end
			SHIFT[t] = nom/denom # Calculate the amount that requires to be shifted
		end
        @constraint(m, dShift[t=1:T],
            Yvars[:pdv][i,t] <= param[:Pd][i,t] - param[:Pd][i,t] * (1-param[:aslDet][i,t,s]) + param[:aslDet][i,t,s] * SHIFT[t])
        @constraint(m, dShift[t=1:T],
            Yvars[:pdv][i,t] >= param[:Pd][i,t] - param[:Pd][i,t] * (1-param[:aslDet][i,t,s]) + param[:aslDet][i,t,s] * SHIFT[t])
	end


	return m, stoMap
end

function climate_benders_subprob_objective(m::JuMP.Model, Yvars::Dict, scenario::scenarioType, stoMap::Dict; kwargs...)

	data = Dict(kwargs)

    M = config.BIGMDOUBLE

	B = data[:package][:param][:B]
	T = data[:package][:param][:T]

	# Setup the objective function
	@objective(m, Min,  M*sum(Yvars[:etap][i,t] + Yvars[:etan][i,t] for i=1:B, t=1:T))

	return m, stoMap
end
