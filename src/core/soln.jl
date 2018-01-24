"""
    This subroutine union a set of designs by taking the maximum design
"""
function union_design(design_pool::Array; kwargs...)

    options = Dict(kwargs)

    K = length(design_pool)
    B,T = size(design_pool[1].pg)

    design = designType()
    design.pg = zeros(B, T)
    design.h = zeros(B, T)

    for k in 1:K
        for b in 1:B
            for t in 1:T
                design.pg[b,t] = max(design_pool[k].pg[b,t], design.pg[b,t])
                design.h[b,t] = max(design_pool[k].h[b,t], design.h[b,t])
            end
        end
    end

    (haskey(options, :param)) && (design.cost = get_design_cost(design, options[:param])[1])

    # If extra parameters is indicated, it means this union is at higher level which the design has been tested with feasibility
    if haskey(options, :S)
        design.feamap = zeros(Int, extras[:S])
        for s in 1:extras[:S]
            for k in 1:K
                # If any one is feasible, then the join is feasible
                design.feamap[s] = max(design_pool[k].feamap[s], 0)
            end
        end
    end

    return design
end

# This subroutine parse the design information from a optimization model
function get_design(model::JuMP.Model)

    design = designType()
    varPg = getindex(model, :pg)
    varH = getindex(model, :h)
    valPg = round.(getvalue(varPg))
    valH = round.(getvalue(varH))
    design.pg = sparse(valPg)
    design.h = sparse(valH)
    design.cost = getobjectivevalue(model)

    return design
end

function get_design_cost(design, param::Dict)

    if isa(design, designType)
        pg = design.pg
        h = design.h
    elseif isa(design, solnType)
        pg = design.primal[:pg]
        h = design.primal[:h]
    else
        error("ERROR|soln.jl|get_design_cost()|Unknown solution type...")
    end

    cost = 0
    expandcost = 0
    hardencost = 0
    B = param[:B]
    T = param[:T]

    # This cost function is problem specific
    for b in 1:B
        expandcost += (pg[b,1] - param[:Pg0][b])*param[:Cg][b,1]
        hardencost += (h[b,1] - param[:H0][b])*param[:Ch][b,1]
        for t in 2:T
            expandcost += (pg[b,t] - pg[b,t-1])*param[:Cg][b,t]
            hardencost += (h[b,t] - h[b,t-1])*param[:Ch][b,t]
        end
    end

    return expandcost + hardencost, expandcost, hardencost
end

function get_master_selection(master::oneProblem, exargs::Dict)

    S = exargs[:S]

    # Find out what scenarios did the master problem picked
    pickScenarioPool = []
    neglectedScenarioPool = []

    # Scan for scenarios that optimal master decision made feasible
    varY = getindex(master.model, :Y)
    varACT = getindex(master.model, :ACT)
    Y = getvalue(varY)
    ACT = getvalue(varACT)

    print("INFO: Master Selection >> @M[")
    for s in 1:S
        if sum(ACT[s]) > 0
            push!(pickScenarioPool, s)
            print("+")
        else
            push!(neglectedScenarioPool, s)
            print("-")
        end
    end
    print("]M@\n")

    return pickScenarioPool, neglectedScenarioPool

end

# This subroutine generate a null design
function null_design(param::Dict)
    null = designType()
    null.pg = zeros(Int, param[:B], param[:T])
    null.h = zeros(Int, param[:B], param[:T])
    return null
end

function infea_design(source, param::Dict)
    infea = designType()
    infea.pg = zeros(Int, param[:B], param[:T])
    infea.h = zeros(Int, param[:B], param[:T])
    infea.active = false
    infea.cost = Inf
    infea.lb = -Inf
    infea.time = 0.0
    infea.source = [source]
    infea.feamap = zeros(Int, param[:S])
    return infea
end


# This subroutine prints out the design in a user-friendly way
function print_design(design, param::Dict)

    allSolution = false

    if isa(design, designType)
        pg = design.pg
        h  = design.h
    elseif isa(design, solnType)
        pg = design.primal[:pg]
        h  = design.primal[:h]
        f  = design.primal[:f]
    else
        errro("ERROR|soln.jl|print_design()|unkown solution type.")
    end

    B, T = size(pg)

    # Collect the sets of bus
    buildBus = Set()
    hardenBus = Set()

    for b in 1:B
        if (pg[b,T] - param[:Pg0][b]) > 0.1
            push!(buildBus,b)
        end
        if (h[b,T] - param[:H0][b]) > 0.1
            push!(hardenBus,b)
        end
    end

    totalExpand = sum(pg[:,T]) - sum(param[:Pg0])
    totalExpandBus = length(buildBus)
    info("Expansion Design [New Generator Added $totalExpand on $totalExpandBus buses]")
    for b in buildBus
        print(string("[EXPAND]BUS $(b)[Ele=$(round.(param[:Ele][b],2)), Cost=$(round.(param[:Cg][b,1],2)), Cap=$(round.(param[:PgUB][b],2)), Top=$(param[:Pgbar][b]), Init=$(param[:Pg0][b])]"))
        for t in 1:T
            print(round.(abs(pg[b,t]), 0)," -> ")
        end
        print("e\n")
    end

    totalHarden = sum(h[:,T]) - Int(sum(param[:H0]))
    totalHardenBus = length(hardenBus)
    info("Harden Design [New harden added $totalHarden on $totalHardenBus buses]")
    for b in hardenBus
        print(string("[HARDEN]BUS $(b)[Ele=$(round.(param[:Ele][b],2)), Cost=$(round.(param[:Ch][b,1],2)), Cap=$(round.(param[:PgUB][b],2)), Top=$(param[:Pgbar][b]), Init=$(param[:Pg0][b])] "))
        for t in 1:T
            print(round.(abs(h[b,t]), 0)," -> ")
        end
        print("e\n")
    end

end

function print_demand_summary(prob::oneProblem, sol::solnType)

    S = prob.param[:S]
    T = prob.param[:T]
    B = prob.param[:B]

    origTimeTotalLoad = []

    for s in 1:S
        for t in 1:T
            timeTotalLoad = round.(sum(prob.param[:aslDetPd][:,t,s]),2)
            info("Total Load at Time Period $t >> $timeTotalLoad ;")
            push!(origTimeTotalLoad, timeTotalLoad)
        end
        for t in 1:T
            timeTotalLoad = round.(sum(sol.primal[:pdv][:,t,s]),2)
            info("[SCEN$s] Total load at time period $t >> $timeTotalLoad [", round.(100*timeTotalLoad/origTimeTotalLoad[t],2) ,"]")
            allShedding = 0.0
            allShifted = 0.0
            origAssLoad = 0.0
            for i in 1:B
                if prob.param[:aslDetPd][i,t,s] - sol.primal[:pdv][i,t,s] > 0.1
                    sheddingLoad = prob.param[:aslDetPd][i,t,s]-sol.primal[:pdv][i,t,s]
                    info("[SCEN$s][Bus$i][--] Losing Load ", sheddingLoad,"]","[",sol.primal[:pdv][i,t,s],"];")
                    @assert abs(sol.primal[:ass][i,t,s] - 0.0) <= 0.1
                    allShedding += sheddingLoad
                elseif prob.param[:aslDetPd][i,t,s] - sol.primal[:pdv][i,t,s] < -0.1
                    shiftedLoad = sol.primal[:pdv][i,t,s]-prob.param[:aslDetPd][i,t,s]
                    info("[SCEN$s][Bus$i][++] Adding Load ", shiftedLoad,"]","[",sol.primal[:pdv][i,t,s],"];")
                    @assert abs(sol.primal[:ass][i,t,s] - 0.0) <= 0.1
                    allShifted += shiftedLoad
                end
                origAssLoad += (1-sol.primal[:ass][i,t,s])*prob.param[:aslDetPd][i,t,s]
            end
            info("Total Shedding Load >> $allShedding; Total Shifted Load $allShifted;")
            info("Shedding Percentage >> ",round.(100*(allShedding-allShifted)/origAssLoad) ,"%;")
        end
    end

end

"""
    This subroutine parse a JSON file into a design type.
"""
function parse_design(dPath::AbstractString, param::Dict)

    df = JSON.parsefile(dPath)

    B = length(df["pg"][1])
    T = length(df["pg"])

    design = designType()
    design.pg = zeros(Int,B,T)
    design.h = zeros(Int,B,T)
    design.feamap = []
    design.cost = 0.0

	if !isa(df["pg"][1][1], Float64)
		return :Infeasible
	end
    for t in 1:T
        for b in 1:B
            design.pg[b,t] = Int(round.(df["pg"][t][b]))
            design.h[b,t] = Int(round.(df["h"][t][b]))
        end
    end

    design.cost, non, non = get_design_cost(design, param)
    return design
end

#This function works with the solutions and construct understandable
#   interpretations of given solutions in Model().
function get_primal_solution(prob::oneProblem, sol::solnType=solnType())

	sol.primal = Dict();

	if isa(prob, oneProblem)
		# Record all primal solution from the model
		for i in keys(prob.vars)
			varI = getindex(prob.model,i)
			sol.primal[i] = getvalue(varI)
		end
	elseif isa(prob, JuMP.Model)
		# Prob itself is an JuMP model, no need to worry
		for i in keys(prob.varDict)
			varI = getindex(prob,i)
			sol.primal[i] = getvalue(prob.varDict[i])
		end
	end

	return sol
end

function warm_start_problem(prob::oneProblem, sol::solnType=solnType())
    for var in keys(oneProblem.vars)
        setvalue(oneProblem.vars[var], sol.primal[var])
    end
    return oneProblem
end


"""
    Currently utilized in deterministic algorithm.
"""
function print_scenario_selection(param::Dict, sol::solnType, exargs::Dict)

    if "sample-based-risk" in exargs[:FEATURES]
        S = param[:S]
        feaPool = []
        for s in 1:S
            if sol.primal[:fs][s] > 0.0
                push!(feaPool, s)
            end
        end
        print("Final Solution Feasible Scenarios $feaPool \n")
    elseif "timm-based-risk" in exargs[:FEATURES]
        T = param[:T]
        S = param[:S]
        println("No implementation yet.")
    end

end

"""
    TODO :
"""
function summary_design(param::Dict, design)

    summary = Dict()

    # TODO: collection information

    return summary

end


"""
    Handle two designs and compare them with details
"""
function compare_designs(param::Dict, solA, solB)

    if isa(solA, designType)
        buildA = solA.pg
        expandA = sol.h
    elseif isa(solA, solnType)
        buildA = solA.primal[:pg]
        expandA = solA.primal[:h]
    else
        error("ERROR|soln.jl|compare_designs|Unkown solution type")
    end

    if isa(solB, designType)
        buildB = solB.pg
        expandB = solB.h
    elseif isa(solA, solnType)
        buildB = solB.primal[:pg]
        expandB = solB.primal[:h]
    else
        error("ERROR|soln.jl|compare_designs|Unkown solution type")
    end

    # Measure the matrix for parameter fetch, all other informaiton stored in param
    B, T = size(buildA)

    # Build comparsion
    #   Total Build
    #   Bus Selected
    #   Location difference
    #   ...

    # Expand comparison
    #   Total Expand
    #   Bus Selected
    #   Location differentce
    #   ...

end
