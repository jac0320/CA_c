function base_formulation(prob::Dict, param::Dict, stoc::stocType, exargs::Dict, selection=[]; kwargs...)

	base = oneProblem()
	param = check_parameter_intactness(param)

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	post_adaptation_vars(base, param)
	post_logical_vars(base, param, selection)

	post_incremental_cons(base, param)
	post_logical_cons(base, param, stoc, selection)
	post_risk_cons(base, param, exargs, selection)

	post_adaptation_obj(base, param)

	return base
end

function cb_model(prob::oneProblem, param::Dict, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)
	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	post_cb_vars(prob, param, selection)
	post_cb_cons(prob, param, selection, subprobType)

	return prob
end

function cnf_model(prob::oneProblem, param::Dict, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)
	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	haskey(options, :logical) ? logical = options[:logical] : logical = true
	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	post_cnf_vars(prob, param, selection, subprobType)
	post_cnf_cons(prob, param, selection, subprobType)

	return prob
end

function dcpf_model(prob::oneProblem, param::Dict, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)

	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	post_dcpf_vars(prob, param, selection, subprobType)
	post_dcpf_cons(prob, param, selection, subprobType)

	return prob
end

function post_adaptation_vars(p::oneProblem, param::Dict)

    B, T = param[:B], param[:T]

    # [TODO] Make var type flexiable
    p.vars[:pg] = @variable(p.model, 0<=pg[i=1:B, 1:T]<=param[:Pgbar][i], Int)
    p.vars[:h] = @variable(p.model, 0<=h[i=1:B, 1:T]<=param[:Hbar][i], Int)

    return
end

function post_logical_vars(p::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T = length(selection), param[:B], param[:T]

    p.vars[:fs] = @variable(p.model, fs[1:S], Bin)  # Added from behind
    p.vars[:f] = @variable(p.model, f[1:T, 1:S], Bin)
    p.vars[:a] = @variable(p.model, a[1:B, 1:T, 1:S], Bin)
    p.vars[:ass] = @variable(p.model, ass[1:B, 1:T, 1:S], Bin)
    p.vars[:ap] = @variable(p.model, ap[1:B, 1:T, 1:S], Int)
    p.vars[:ah] = @variable(p.model, ah[1:B, 1:T, 1:S], Int)

    return
end

function post_incremental_cons(p::oneProblem, param::Dict)

    B, T = param[:B], param[:T]

    @constraint(p.model, [i=1:B], p.vars[:pg][i,1] >= param[:Pg0][i])
    @constraint(p.model, [i=1:B], p.vars[:h][i,1] >= param[:H0][i])
    @constraint(p.model, [i=1:B,t=2:T], p.vars[:pg][i,t-1] <= p.vars[:pg][i,t])
    @constraint(p.model, [i=1:B,t=2:T], p.vars[:h][i,t-1] <= p.vars[:h][i,t])

    return
end

function post_logical_cons(prob::oneProblem, param::Dict, sto::stocType, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T = length(selection), param[:B], param[:T]

    # Above are reduced due to the information can be determinied deterministically
    @constraint(prob.model, [i=1:B, t=1:T, s=1:S],
        sto.scenarios[s].data["SS"][i,t] * (2*prob.vars[:ass][i,t,s]-1) -
            2*param[:Ele][i]*prob.vars[:ass][i,t,s] -
            2*param[:ProM][i]*prob.vars[:ah][i,t,s] +
            param[:Ele][i] + param[:ProM][i]*prob.vars[:h][i,t] <= 0)

    # McCormick Relaxation => ap = a * pg # Replace this moduleswith methods in relax.jl
    for i in 1:B
        for t in 1:T
            for s in 1:S
                @constraint(prob.model, prob.vars[:a][i,t,s] == param[:aslDet][i,t,selection[s]] * prob.vars[:ass][i,t,s])
                mccormick(prob.model, prob.vars[:ap][i,t,s], prob.vars[:a][i,t,s], prob.vars[:pg][i,t], 0, 1, param[:Pg0][i], param[:Pgbar][i]);
                mccormick(prob.model, prob.vars[:ah][i,t,s] ,prob.vars[:ass][i,t,s], prob.vars[:h][i,t], 0, 1, param[:H0][i], param[:Hbar][i]);
            end
        end
    end

    return
end

function post_adaptation_obj(p::oneProblem, param::Dict)

    B, T = param[:B], param[:T]

    @objective(p.model, Min, sum(param[:Cg][i,1]*(p.vars[:pg][i,1]-param[:Pg0][i]) for i=1:B)
                        + sum(param[:Ch][i,1]*(p.vars[:h][i,1]-param[:H0][i]) for i=1:B)
                        + sum(param[:Cg][i,t]*(p.vars[:pg][i,t]-p.vars[:pg][i,t-1]) for i=1:B, t=2:T)
                        + sum(param[:Ch][i,t]*(p.vars[:h][i,t]-p.vars[:h][i,t-1]) for i=1:B, t=2:T));

    return
end

function post_cb_vars(p::oneProblem, param::Dict, selection=[], sbtype="free")

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, B, T = length(selection), param[:B], param[:T]

    @assert length(p.vars[:fs]) == S

    p.vars[:pdv] = @variable(p.model, pdv[1:B, 1:T, 1:S]>=0)

    if sbtype == "slackness"
        p.vars[:scapslack] = @variable(p.model, capslack[1:T,1:S]>=0)
    end

    return
end

function post_cb_cons(p::oneProblem, param::Dict, selection=[], sbtype="free")


    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, B, T = length(selection), param[:B], param[:T]

    if sbtype == "free"
        @constraint(p.model, [t=1:T,s=1:S],
            sum(param[:PgUB][i] * p.vars[:ap][i,t,s] for i=1:B) -
            sum(p.vars[:pdv][i,t,s] for i=1:B) >=
            -sum(param[:aslDetPd][i,t,s] for i=1:B)*(1-p.vars[:f][t,s]))
        @constraint(p.model, [t=1:T,s=1:S],
            sum(param[:aslDetPd][i,t,s] for i=1:B) * p.vars[:f][t,s] >=
            sum(param[:PgUB][i] * p.vars[:ap][i,t,s] for i=1:B) -
            sum(p.vars[:pdv][i,t,s] for i=1:B))
    elseif sbtype == "tight"
        @constraint(p.model, [t=1:T,s=1:S],
            sum(param[:PgUB][i] * p.vars[:ap][i,t,s] for i=1:B) >= sum(p.vars[:pdv][i,t,s] for i=1:B))
    elseif sbtype == "slackness"
        @assert haskey(p.vars, :capslack)
        @constraint(p.model, [t=1:T,s=1:S],
            sum(param[:PgUB][i] * p.vars[:ap][i,t,s] for i=1:B) +
            p.vars[:capslack][t,s] >= sum(param[:Pd][i,t] for i=1:B))
        post_cb_slack_obj(p)
    else
        error("ERROR|general.jl|cb_model()|Unknown subproblem type.")
    end

    post_cbloadshed_cons(p, param, selection)

    return
end

function post_cb_slack_obj(prob::oneProblem)
    @objective(prob.model, Min, sum(prob.vars[:capslack]))
    return
end

function post_cnf_slack_obj(prob::oneProblem)

    @objective(prob.model, Min,
        sum(prob.vars[:slsSlack]) + sum(prob.vars[:flowSlackpos]) + sum(prob.vars[:flowSlackneg]))

    return
end

function post_dcpf_slack_obj(prob::oneProblem)

    @objective(prob.model, Min, sum(prob.vars[:slsSlack]) +
                                sum(prob.vars[:flowSlackpos]) +
                                sum(prob.vars[:flowSlackneg]) +
                                sum(prob.vars[:flowEquSlackpos]) +
                                sum(prob.vars[:flowEquSlackneg]))

    return
end

function post_cnf_vars(prob::oneProblem, param::Dict, selection=[], sbtype="free")

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    #Power Flow i->j | Actual Generation | Flow Balance Feasible Indicator

    # If a bus is flooded, then it linked thermal limits will be cutted off
    # If a bus is undersea permentaly, also cut off the associated thermal limits
    prob.vars[:p] = @variable(prob.model, p[i=1:B,j=1:B,1:T,1:S; param[:EDGE][i,j]==1])
    prob.vars[:pdv] = @variable(prob.model, pdv[1:B, 1:T, 1:S]>=0) # Variable Demand
    prob.vars[:psg] = @variable(prob.model, psg[1:B, 1:T, 1:S]>=0) # Variable generation output
    prob.vars[:fn] = @variable(prob.model, fn[1:B, 1:T, 1:S], Bin)

    if sbtype == "slackness"
        prob.vars[:flowSlackpos] = @variable(prob.model, flowSlackpos[1:B,1:T,1:S]>=0)
        prob.vars[:flowSlackneg] = @variable(prob.model, flowSlackneg[1:B,1:T,1:S]>=0)
    end

    return
end

function post_switching_cons(prob::oneProblem, param::Dict, selection=[])
    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, B, T = length(selection), param[:B],param[:T]

    for i in 1:B
        # Line switch off given sea level rise
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["in"]],
            prob.vars[:p][j,i,t,s] <= prob.vars[:ass][i,t,s] * param[:Lcap][j,i])
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["in"]],
            prob.vars[:p][j,i,t,s] >= prob.vars[:ass][i,t,s] * (-param[:Lcap][j,i]))
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["out"]],
            prob.vars[:p][i,j,t,s] <= prob.vars[:ass][i,t,s] * param[:Lcap][i,j])
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["out"]],
            prob.vars[:p][i,j,t,s] >= prob.vars[:ass][i,t,s] * (-param[:Lcap][i,j]))
        # Line swtich off given storm surge
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["in"]],
            prob.vars[:p][j,i,t,s] <= param[:aslDet][i,t,selection[s]] * param[:Lcap][j,i])
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["in"]],
            prob.vars[:p][j,i,t,s] >= param[:aslDet][i,t,selection[s]] * (-param[:Lcap][j,i]))
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["out"]],
            prob.vars[:p][i,j,t,s] <= param[:aslDet][i,t,selection[s]] * param[:Lcap][i,j])
        @constraint(prob.model, [j=1:B, t=1:T, s=1:S; j in param[:Edge][i]["out"]],
            prob.vars[:p][i,j,t,s] >= param[:aslDet][i,t,selection[s]] * (-param[:Lcap][i,j]))
    end

    return
end

function post_generation_cons(prob::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    # Varibale generation ouput upper bounded by total capaciated units
    @constraint(prob.model, [i=1:B,t=1:T,s=1:S],
        prob.vars[:psg][i,t,s] <= param[:PgUB][i] * prob.vars[:ap][i,t,s])

    return
end

function post_cnf_logic_cons(prob::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    @constraint(prob.model, [s=1:S,t=1:T],
        B*(prob.vars[:f][t,s]-1) <= sum(prob.vars[:fn][i,t,s] for i=1:B))
    @constraint(prob.model, [s=1:S,t=1:T],
        sum(prob.vars[:fn][i,t,s] for i=1:B) <= B*(1-prob.vars[:f][t,s]))

    return
end

function post_flowbalance_cons(prob::oneProblem, param::Dict, selection=[], sbtype="free", post_additional_logic_cons::Function=nothing)

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    if sbtype == "free"
        @assert haskey(prob.vars, :f)
        for i in 1:B
            @constraint(prob.model, [t=1:T,s=1:S],
                sum(prob.vars[:p][j,i,t,s] for j in param[:Edge][i]["in"]) -
                sum(prob.vars[:p][i,j,t,s] for j in param[:Edge][i]["out"]) +
                prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s] <= 999999 * prob.vars[:fn][i,t,s])
            @constraint(prob.model, [t=1:T,s=1:S],
                -(sum(prob.vars[:p][j,i,t,s] for j in param[:Edge][i]["in"]) -
                sum(prob.vars[:p][i,j,t,s] for j in param[:Edge][i]["out"]) +
                prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s]) <= 999999 * prob.vars[:fn][i,t,s])
        end
        post_additional_logic_cons(prob, param, selection)
    elseif sbtype == "tight"
        for i in 1:B
            @constraint(prob.model, [t=1:T, s=1:S],
                sum(prob.vars[:p][j,i,t,s] for j in param[:Edge][i]["in"]) - sum(prob.vars[:p][i,j,t,s] for j in param[:Edge][i]["out"]) + prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s] == 0)
        end
    elseif sbtype == "slackness"
        for i in 1:B
            @constraint(prob.model, [t=1:T,s=1:S],
                sum(prob.vars[:p][j,i,t,s] for j in param[:Edge][i]["in"]) - sum(prob.vars[:p][i,j,t,s] for j in param[:Edge][i]["out"]) + prob.vars[:psg][i,t,s] - prob.vars[:pdv][i,t,s] + prob.vars[:flowSlackpos][i,t,s] - prob.vars[:flowSlackneg][i,t,s] == 0)
        end
        post_cnf_slack_obj(prob)
    else
        error("ERROR|general.jl|cnf_model()|Unkown subproblem type.")
    end
end

function post_cnf_cons(prob::oneProblem, param::Dict, selection=[], sbtype="free")

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    post_switching_cons(prob, param, selection)
    post_generation_cons(prob, param, selection)
    post_flowbalance_cons(prob, param, selection, sbtype, post_cnf_logic_cons)
    post_loadshed_cons(prob, param, selection)

    return
end

function post_dcpf_vars(prob::oneProblem, param::Dict, selection=[], sbtype="free")

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    post_cnf_vars(prob, param, selection, sbtype)
    prob.vars[:theta] = @variable(prob.model,param[:AngleLimit]>=theta[1:B, 1:T, 1:S]>=-param[:AngleLimit])
    prob.vars[:assss] = @variable(prob.model, assss[i=1:B,j=1:B,1:T,1:S; param[:EDGE][i,j] == 1], Bin)

    if sbtype == "free"
        prob.vars[:fln] = @variable(prob.model, fln[i=1:B,j=1:B,1:T,1:S; param[:EDGE][i,j] == 1], Bin)
    elseif sbtype == "slackness"
        prob.vars[:flowEquSlackpos] = @variable(prob.model,
            flowEquSlackpos[i=1:B,j=1:B,t=1:T,s=1:S; param[:EDGE][i,j] == 1]>=0)
        prob.vars[:flowEquSlackneg] = @variable(prob.model,
            flowEquSlackneg[i=1:B,j=1:B,t=1:T,s=1:S; param[:EDGE][i,j] == 1]>=0)
    else
        error("")
    end

    return
end

function post_dcpf_cons(prob::oneProblem, param::Dict, selection=[], sbtype="free")

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    post_switching_cons(prob, param, selection)
    post_generation_cons(prob, param, selection)
    post_flowbalance_cons(prob, param, selection, sbtype, post_dcpf_logic_cons)

    post_refbus_cons(prob, param, selection)
    post_anglelimits_cons(prob, param, selection)
    post_assss_cons(prob, param, selection)
    post_flowequation_cons(prob, param, selection, sbtype)
    post_loadshed_cons(prob, param, selection)

    return
end

function post_dcpf_logic_cons(prob::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    # Add Logical Feasibility: this means the logical constraints in network functin is passed
    @constraint(prob.model, [s=1:S,t=1:T],
        (B+L)*(prob.vars[:f][t,s]-1) <= sum(prob.vars[:fn][i,t,s] for i=1:B) + sum(prob.vars[:fln][i,j,t,s] for i=1:B,j=1:B if param[:EDGE][i,j]==1))
    @constraint(prob.model, [s=1:S,t=1:T],
        sum(prob.vars[:fn][i,t,s] for i=1:B) + sum(prob.vars[:fln][i,j,t,s] for i=1:B, j=1:B if param[:EDGE][i,j]==1) <= (B+L)*(1-prob.vars[:f][t,s]))

    return
end

function post_anglelimits_cons(prob::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    @constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; param[:EDGE][i,j] == 1],
        prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] <= param[:AngleShiftLimit])
    @constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; param[:EDGE][i,j] == 1],
        prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s] >= -param[:AngleShiftLimit])

    return
end

function post_assss_cons(prob::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    # A huge extra system of constraints if DC model is included
    # [TODO] write better loop
    for i in 1:B
        for j in 1:B
            if param[:EDGE][i,j] == 1
                for t in 1:T
                    for s in 1:S
                        mccormick(prob.model, prob.vars[:assss][i,j,t,s], prob.vars[:ass][i,t,s], prob.vars[:ass][j,t,s], 0, 1, 0, 1)
                    end
                end
            end
        end
    end

    return
end

function post_refbus_cons(prob::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    @constraint(prob.model, [s=1:S, t=1:T], prob.vars[:theta][param[:RefBus],t,s] == 0)

    return
end

function post_flowequation_cons(prob::oneProblem, param::Dict, selection=[], sbtype="free")

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    if sbtype == "free"
        @constraint(prob.model, [t=1:T,s=1:S, i=1:B, j=1:B; param[:EDGE][i,j] == 1], # Seal Upper bound
            prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/param[:lineX][i,j] <=
            (param[:Lcap][i,j] + param[:AngleShiftLimit] / param[:lineX][i,j]) *
            prob.vars[:fln][i,j,t,s] + 100*(param[:Lcap][i,j] + param[:AngleShiftLimit] /
            param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] *
            param[:aslDet][i,t,selection[s]] * param[:aslDet][j,t,selection[s]]))
        @constraint(prob.model, [t=1:T,s=1:S, i=1:B, j=1:B; param[:EDGE][i,j] == 1], # Seal lower bound
            -(prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] - prob.vars[:theta][j,t,s])/param[:lineX][i,j]) <=
              (param[:Lcap][i,j] + param[:AngleShiftLimit] / param[:lineX][i,j]) *
              prob.vars[:fln][i,j,t,s] + 100*(param[:Lcap][i,j] + param[:AngleShiftLimit] /
              param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] *
              param[:aslDet][i,t,selection[s]] * param[:aslDet][j,t,selection[s]]))

    elseif sbtype == "tight"
        @constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; param[:EDGE][i,j] == 1], # Seal Upper bound
            prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] -
            prob.vars[:theta][j,t,s])/param[:lineX][i,j] <=
            100*(param[:Lcap][i,j] + param[:AngleShiftLimit] /
            param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] *
            param[:aslDet][i,t,selection[s]] * param[:aslDet][j,t,selection[s]]))
        @constraint(prob.model, [t=1:T, s=1:S, i=1:B, j=1:B; param[:EDGE][i,j] == 1], # Seal lower bound
           -(prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] -
           prob.vars[:theta][j,t,s])/param[:lineX][i,j]) <=
           100*(param[:Lcap][i,j] + param[:AngleShiftLimit] /
           param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] *
           param[:aslDet][i,t,selection[s]] * param[:aslDet][j,t,selection[s]]))

    elseif sbtype == "slackness"
        @constraint(prob.model, [t=1:T,s=1:S,i=1:B,j=1:B; param[:EDGE][i,j]==1], # Seal Upper bound
            prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] -
            prob.vars[:theta][j,t,s])/param[:lineX][i,j] +
            prob.vars[:flowEquSlackpos][i,j,t,s] - prob.vars[:flowEquSlackneg][i,j,t,s] <=
            (param[:Lcap][i,j] + param[:AngleShiftLimit] /
            param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] *
            param[:aslDet][i,t,selection[s]] * param[:aslDet][j,t,selection[s]]) )
        @constraint(prob.model, [t=1:T,s=1:S,i=1:B,j=1:B; param[:EDGE][i,j]==1], # Seal lower bound
            -(prob.vars[:p][i,j,t,s] - (prob.vars[:theta][i,t,s] -
            prob.vars[:theta][j,t,s])/param[:lineX][i,j] +
            prob.vars[:flowEquSlackpos][i,j,t,s] - prob.vars[:flowEquSlackneg][i,j,t,s]) <=
            (param[:Lcap][i,j] + param[:AngleShiftLimit] /
            param[:lineX][i,j]) * (1 - prob.vars[:assss][i,j,t,s] *
            param[:aslDet][i,t,selection[s]] * param[:aslDet][j,t,selection[s]]) )
        post_dcpf_slack_obj(prob)

    else
        error("ERROR|general.jl|dcpf_model()|Unkown subprob type.")
    end

    return
end

function post_risk_cons(prob::oneProblem, param::Dict, driver::Dict, selection=[]; overrideEps= -1.0)

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]
    overrideEps > 0.0 ? eps = overrideEps : eps = driver[:eps]

    @constraint(prob.model, [s=1:S],
        T * prob.vars[:fs][s] <= sum(prob.vars[:f][t,s] for t=1:T))
    @constraint(prob.model, risk,
        sum(prob.vars[:fs][s] for s=1:S) >= ceil(S * (1 - eps)))

    return
end

function post_cbloadshed_cons(prob::oneProblem, param::Dict, selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    @constraint(prob.model, [t=1:T, s=1:S],
        sum(prob.vars[:pdv][i,t,s] for i=1:B) >= param[:SHEDLambda] *
        sum(param[:aslDetPd][i,t,selection[s]] for i=1:B))
    @constraint(prob.model, dispLb[i=1:B, t=1:T, s=1:S],
        prob.vars[:pdv][i,t,s] <= param[:aslDetPd][i,t,selection[s]]*prob.vars[:ass][i,t,s])

    return
end

function post_loadshed_cons(prob::oneProblem, param::Dict, selection=[], sbtype="free")

    isempty(selection) ? selection = [1:param[:S];] : selection = selection
    S, B, T, L = length(selection), param[:B], param[:T], param[:L]

    if sbtype == "free"
        @assert haskey(prob.vars, :f)
        @constraint(prob.model, [t=1:T, s=1:S],
            sum(prob.vars[:pdv][i,t,s] for i=1:B) >= param[:SHEDLambda] *
            sum(param[:aslDetPd][i,t,selection[s]] for i=1:B) * prob.vars[:fs][s])
        @constraint(prob.model, dispLb[i=1:B, t=1:T, s=1:S],
            prob.vars[:pdv][i,t,s] <= param[:aslDetPd][i,t,selection[s]] * prob.vars[:ass][i,t,s])

    elseif sbtype == "slackness"
        prob.vars[:slsSlack] = @variable(prob.model, slsslack[1:T,1:S]>=0)
        @constraint(prob.model, shed[t=1:T, s=1:S],
            sum(prob.vars[:pdv][i,t,s] for i=1:B) + prob.vars[:slsSlack][t,s] >=
            param[:SHEDLambda] * sum(param[:aslDetPd][i,t,selection[s]] for i=1:B))
        @constraint(prob.model, dispLb[i=1:B, t=1:T, s=1:S],
            prob.vars[:pdv][i,t,s] <= param[:aslDetPd][i,t,selection[s]] * prob.vars[:ass][i,t,s])

    else sbtype == "tight"
        @constraint(prob.model, [i=1:B, t=1:T, s=1:S],
            prob.vars[:pdv][i,t,s] <= param[:aslDetPd][i,t,selection[s]] * prob.vars[:ass][i,t,s])
        @constraint(prob.model, [t=1:T, s=1:S],
            sum(prob.vars[:pdv][i,t,s] for i=1:B) >= param[:SHEDLambda] *
            sum(param[:aslDetPd][i,t,selection[s]] for i=1:B))
    end

    return
end

function mccormick(m,xy,x,y,xˡ,xᵘ,yˡ,yᵘ)
    @constraint(m, xy >= xˡ*y + yˡ*x - xˡ*yˡ)
    @constraint(m, xy >= xᵘ*y + yᵘ*x - xᵘ*yᵘ)
    @constraint(m, xy <= xˡ*y + yᵘ*x - xˡ*yᵘ)
    @constraint(m, xy <= xᵘ*y + yˡ*x - xᵘ*yˡ)
    return m
end
