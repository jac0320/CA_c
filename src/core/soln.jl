function union_design(pool::Array, param::Dict; source::Array=[])

    P, B, T, S = length(pool), param[:B], param[:T], param[:S]

    pg = zeros(Int, B, T)
    h = zeros(Int, B, T)
    f = zeros(Int, param[:S])

    for k in 1:P
        for b in 1:B, t in 1:T
            pg[b,t] = max(pool[k].pg[b,t], pg[b,t])
            h[b,t] = max(pool[k].h[b,t], h[b,t])
        end
        for s in 1:S
            f[s] = max(pool[k].feamap[s], f[s])
        end
    end

    id = length(pool) + 1

    d = designType(id, pg, h)
    d.feamap = f
    d.cost = get_design_cost(pg, h, param)
    d.coverage = sum(d.feamap) / length(d.feamap)
    d.source = source

    return d
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

function get_design(prob::oneProblem)

    valPg = convert(Array{Int}, round.(getvalue(prob.vars[:pg])))
    valH = convert(Array{Int}, round.(getvalue(prob.vars[:h])))

    d = designType(idx, sparse(valPg), sparse(valH), getobjectivevalue(prob.model))
    d.lb = prob.model.objBound
    d.time = getsolvetime(prob.model)

    return d
end

function get_design_cost(pg::Array, h::Array, param::Dict)

    expandcost = 0
    hardencost = 0
    B, T = param[:B], param[:T]

    # This cost function is problem specific
    for b in 1:B
        expandcost += (pg[b,1] - param[:Pg0][b])*param[:Cg][b,1]
        hardencost += (h[b,1] - param[:H0][b])*param[:Ch][b,1]
        for t in 2:T
            expandcost += (pg[b,t] - pg[b,t-1])*param[:Cg][b,t]
            hardencost += (h[b,t] - h[b,t-1])*param[:Ch][b,t]
        end
    end

    return expandcost + hardencost
end

function get_design_cost(design::Any, param::Dict)

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

function infea_design(source::Any, param::Dict)
    infea = designType()
    infea.pg = zeros(Int, param[:B], param[:T])
    infea.h = zeros(Int, param[:B], param[:T])
    infea.active = false
    infea.cost = Inf
    infea.lb = -Inf
    infea.time = 0.0
    infea.source = collect(source)
    infea.feamap = zeros(Int, param[:S])
    return infea
end

function collect_design(pool::Vector{designType}, d::designType; id=-1)

    id > 0 ? d.k = id : d.k = d.k
    show_design(d)
    push!(pool, d)

    return
end



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

    totalcost, expandcost, hardencost = get_design_cost(design, param)
    info("[COST] $totalcost = $expandcost + $hardencost")

    B, T = size(pg)
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
    info("Expansion || UNIT $totalExpand || BUS $totalExpandBus")
    for b in buildBus
        print("[EXPAND]BUS $(b)[Ele=$(round.(param[:Ele][b],2)), Cost=$(round.(param[:Cg][b,1],2)), Cap=$(round.(param[:PgUB][b],2)), Top=$(param[:Pgbar][b]), Init=$(param[:Pg0][b])]")
        for t in 1:T
            print(round.(abs(pg[b,t]), 0)," -> ")
        end
        print("e\n")
    end

    totalHarden = sum(h[:,T]) - Int(sum(param[:H0]))
    totalHardenBus = length(hardenBus)
    info("HARDENING || UNIT $totalHarden || BUS $totalHardenBus")
    for b in hardenBus
        print("[HARDEN]BUS $(b)[Ele=$(round.(param[:Ele][b],2)), Cost=$(round.(param[:Ch][b,1],2)), Cap=$(round.(param[:PgUB][b],2)), Top=$(param[:Pgbar][b]), Init=$(param[:Pg0][b])] ")
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
