function reactor(power::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

    options = Dict(kwargs)

    T = param[:T]
    B = param[:B]

    # Generate a no hurricane scenario
    info("[REACTOR]Generating null scenario...")
    zeroStoc = null_scenario_stoc(param[:B], param[:T])

    # Solve the no hurricane scenario for a expansion plan
    expandProb = oneProblem()
    expandProb = sbd_base_formulation(power, param, zeroStoc)
    expandProb = attach_scenario(expandProb, zeroStoc, [1], exargs[:MODEL], 0.0, exargs)
    info("[REACTOR]Checking expansion plan using optimization...")
    status = solve(expandProb.model, suppress_warnings=true)

    # Collect results (Dimension [B, T])
    info("[REACTOR]Heuristic expansion plan generated with total cost $(getobjectivevalue(expandProb.model))")
    expandPlan = getvalue(expandProb.vars[:pg])

    # Consider only one scenario of the stocType by random
    pickIdx = rand([1:stoc.S;])
    pickScenario = stoc.scenarios[pickIdx]
    info("[REACTOR]Constructing Hardening base on scenario $(pickIdx)")

    closePlan = ones(B, T)
    for b in 1:B
        for t in 1:T
            if pickScenario.data["SL"][t] >= param[:Ele][b]
                closePlan[b,t] = 0
            end
        end
    end

    hardenPlan = zeros(B, T)

    # Then perform the reactive building
    for t in 2:T
        for b in 1:B
            hardenPlan[b,t] = hardenPlan[b,t-1]
            if closePlan[b,t] == 1 && pickScenario.data["SS"][b,t-1] > param[:Ele][b]
                # Build harden facility incrementally base-on previous time step information (reactively)
                hardenPlan[b,t] = max(hardenPlan[b,t-1], ceil((pickScenario.data["SS"][b,t-1] - param[:Ele][b])/param[:ProM][b]))
                info("[REACTOR]Bus $b need $(hardenPlan[b,t]) hardenings at time step $t")
            end
        end
    end

    # Build hardenning at first step by discount
    for b in 1:B
        hardenPlan[b, 1] = floor(hardenPlan[b, 2]/2)
        if hardenPlan[b, 1] > 0.0
            info("[REACTOR]Bus $b need $(hardenPlan[b,1]) hardenings at time step 1")
        end
    end

    hDesign = designType()
    hDesign.pg = expandPlan
    hDesign.h = hardenPlan

    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info(string("[REACTOR]The total cost is $totalcost = $expandcost + $hardencost"))

    write_output_files(power, param, stoc, hDesign, exargs)
    print_design(hDesign, param)

    return expandProb, hDesign
end

function highland(power::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

    options = Dict(kwargs)

    T = param[:T]
    B = param[:B]

    # Generate a no hurricane scenario
    info("[HIGHLAND] Generating null scenario...")
    zeroStoc = null_scenario_stoc(param[:B], param[:T])

    # Solve the no hurricane scenario for a expansion plan
    expandProb = oneProblem()
    expandProb = sbd_base_formulation(power, param, zeroStoc)
    expandProb = attach_scenario(expandProb, zeroStoc, [1], exargs[:MODEL], 0.0, exargs)

    # meanEle = mean(param[:Ele])
    pickIdx = rand([1:stoc.S;])
    pickScenario = stoc.scenarios[pickIdx]
    percSS = get_percentile_SS(stoc, 90, [pickIdx])

    # percSS = mean(param[:Ele])/2
    info("[HIGHLAND] Criteria 90 percentile recorded storm surage $(percSS)")
    rejectbus = 0
    for i=1:B
        if param[:Ele][i] < percSS
            @constraint(expandProb.model, [t=1:T],
                expandProb.vars[:pg][i,t] <= param[:Pg0][i])
            rejectbus += 1
        end
    end
    info("[HIGHLAND] Rejected total $(rejectbus) buses' location")

    info("[HIGHLAND] Checking expansion plan using optimization...")
    status = solve(expandProb.model, suppress_warnings=true)
    if status == :Infeasible
        info("[HIGHLAND] MARK Infeasible expansion plan located.")
    end

    # Collect results (Dimension [B, T])
    info("[HIGHLAND] Heuristic expansion plan generated with total cost $(getobjectivevalue(expandProb.model))")
    expandPlan = getvalue(expandProb.vars[:pg])

    info("[HIGHLAND] Constructing Hardening base on scenario $(pickIdx)")

    closePlan = ones(B, T)
    for b in 1:B
        for t in 1:T
            if pickScenario.data["SL"][t] >= param[:Ele][b]
                closePlan[b,t] = 0
            end
        end
    end

    hardenPlan = zeros(B, T)
    for t in 2:T
        for b in 1:B
            hardenPlan[b,t] = hardenPlan[b,t-1]
            if closePlan[b,t] == 1 && pickScenario.data["SS"][b,t-1] > param[:Ele][b]
                hardenPlan[b,t] = max(hardenPlan[b,t-1], ceil((pickScenario.data["SS"][b,t-1] - param[:Ele][b])/param[:ProM][b]))
                info("[HIGHLAND] Bus $b need $(hardenPlan[b,t]) hardenings at time step $t")
            end
        end
    end

    for b in 1:B
        hardenPlan[b, 1] = floor(hardenPlan[b, 2]/2)
        if hardenPlan[b, 1] > 0.0
            info("[REACTOR] Bus $b need $(hardenPlan[b,1]) hardenings at time step 1")
        end
    end

    hDesign = designType()
    hDesign.pg = expandPlan
    hDesign.h = hardenPlan

    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info(string("[HIGHLAND] The total cost is $totalcost = $expandcost + $hardencost"))

    write_output_files(power, param, stoc, hDesign, exargs)
    print_design(hDesign, param)

    return expandProb, hDesign
end

function bathtub(power::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

    options = Dict(kwargs)

    T = param[:T]
    B = param[:B]

    # Generate a no hurricane scenario
    info("[BATHTUB]Generating null scenario...")
    zeroStoc = null_scenario_stoc(B, T)

    # Solve the no hurricane scenario for a expansion plan
    expandProb = oneProblem()
    expandProb = sbd_base_formulation(power, param, zeroStoc)
    expandProb = attach_scenario(expandProb, zeroStoc, [1], exargs[:MODEL], 0.0, exargs)

    info("[BATHTUB]Checking expansion plan using basic optimization...")
    status = solve(expandProb.model, suppress_warnings=true)

    # Collect results (Dimension [B, T])
    info("[BATHTUB]Heuristic expansion plan generated with total cost $(getobjectivevalue(expandProb.model))")
    expandPlan = getvalue(expandProb.vars[:pg])

    # Consider only one scenario of the stocType by random
    pickIdx = rand([1:stoc.S;])
    pickScenario = stoc.scenarios[pickIdx]
    maxSS = get_percentile_SS(stoc, 100, [pickIdx])
    info("[BATHTUB]Constructing Hardening base on scenario $(pickIdx)")

    closePlan = ones(B, T)
    for b in 1:B
        for t in 1:T
            if pickScenario.data["SL"][t] >= param[:Ele][b]
                info("[BATHTUB] Bus $b closed due to sea level rise")
                closePlan[b,t] = 0  #Bus is closed due to sea level rise
            end
        end
    end

    hardenPlan = zeros(B, T)
    for t in 1:T
        for b in 1:B
            if t >= 2
                hardenPlan[b,t] = hardenPlan[b,t-1]
            end
            if closePlan[b,t] == 1 && pickScenario.data["SS"][b,t] >= param[:Ele][b]
                # hardenPlan[b,t] = round.((maximum(param[:Ele]) - param[:Ele][b])/param[:ProM][b])
                hardenPlan[b,t] = ceil.((maxSS - param[:Ele][b])/param[:ProM][b])
                info("[BATHTUB]Bus $b need $(hardenPlan[b,t]) hardenings at time step $t")
            end
        end
    end

    hDesign = designType()
    hDesign.pg = expandPlan
    hDesign.h = hardenPlan

    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info(string("[BATHTUB]The total cost is $totalcost = $expandcost + $hardencost"))

    write_output_files(power, param, stoc, hDesign, exargs)
    print_design(hDesign, param)

    return expandProb, hDesign
end

function extreme(power::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

    options = Dict(kwargs)

    T = param[:T]
    B = param[:B]

    # Consider only one scenario of the stocType by random
    pickIdx = rand([1:stoc.S;])
    pickScenario = stoc.scenarios[pickIdx]
    maxSS = get_percentile_SS(stoc, 100, [pickIdx])
    info("[EXTREME]Constructing Hardening base on scenario $(pickIdx)")

    # Generate a no hurricane scenario
    info("[EXTREME]Generating null scenario...")
    zeroStoc = null_scenario_stoc(B, T)
    zeroStoc.scenarios[1].data["SL"] = copy(pickScenario.data["SL"])

    # Solve the no hurricane scenario for a expansion plan
    expandProb = oneProblem()
    expandProb = sbd_base_formulation(power, param, zeroStoc)
    expandProb = attach_scenario(expandProb, zeroStoc, [1], exargs[:MODEL], 0.0, exargs)

    info("[EXTREME]Checking expansion plan using basic optimization...")
    solver_config(expandProb.model, showlog=1, presolve=1, mipgap=0.0)
    status = solve(expandProb.model, suppress_warnings=true)

    writeLP(expandProb.model,"95p.lp")

    # Collect results (Dimension [B, T])
    info("[EXTREME]Heuristic expansion plan generated with total cost $(getobjectivevalue(expandProb.model))")
    expandPlan = getvalue(expandProb.vars[:pg])

    closePlan = ones(B, T)
    for b in 1:B
        for t in 1:T
            if pickScenario.data["SL"][t] >= param[:Ele][b]
                info("[EXTREME] Bus $b closed due to sea level rise")
                closePlan[b,t] = 0  #Bus is closed due to sea level rise
                @show param[:aslDet][b,t]
            else
                @assert param[:aslDet][b,t] == true
            end
        end
    end

    hardenPlan = zeros(B, T)
    for t in 1:T
        for b in 1:B
            if t >=2
                hardenPlan[b,t] = hardenPlan[b,t-1]
            end
            if pickScenario.data["SS"][b,t] >= param[:Ele][b]
                hardenPlan[b,t] = ceil((1.1*maximum(pickScenario.data["SS"]) - param[:Ele][b])/param[:ProM][b])
                # hardenPlan[b,t] = round.(maxSS - param[:Ele][b])/param[:ProM][b]
                info("[EXTREME]Bus $b need $(hardenPlan[b,t]) hardenings at time step $t")
            end
        end
    end

    hDesign = designType()
    hDesign.pg = expandPlan
    hDesign.h = hardenPlan

    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info(string("[EXTREME]The total cost is $totalcost = $expandcost + $hardencost"))

    write_output_files(power, param, stoc, hDesign, exargs)
    print_design(hDesign, param)

    return expandProb, hDesign
end

function warmstart_heuristic(prob::oneProblem, power::Dict, param::Dict, stoc::stocType, exargs::Dict; selection=[])

    T = param[:T]
    B = param[:B]

    isempty(selection) ? selection = [1:stoc.S;] : selection = selection

    # Generate a no hurricane scenario
    zeroStoc = null_scenario_stoc(B, T)

    # Solve the no hurricane scenario for a expansion plan
    expandProb = oneProblem()
    expandProb = sbd_base_formulation(power, param, zeroStoc)
    expandProb = attach_scenario(expandProb, zeroStoc, [1], exargs[:MODEL], 0.0, exargs)

    status = solve(expandProb.model, suppress_warnings=true)

    # Collect results (Dimension [B, T])
    expandPlan = getvalue(expandProb.vars[:pg])

    orig_u_expand = ones(B, T)
    orig_l_expand = ones(B, T)
    orig_u_harden = ones(B, T)
    orig_l_harden = ones(B, T)

    for b in 1:B
        for t in 1:T
            orig_u_expand[b,t] = round.(getupperbound(prob.vars[:pg][b,t]))
            orig_l_expand[b,t] = round.(getlowerbound(prob.vars[:pg][b,t]))
            orig_u_harden[b,t] = round.(getupperbound(prob.vars[:h][b,t]))
            orig_l_harden[b,t] = round.(getlowerbound(prob.vars[:h][b,t]))
            setlowerbound(prob.vars[:pg][b,t], expandPlan[b,t])
            setlowerbound(prob.vars[:h][b,t], orig_u_harden[b,t])
        end
    end

    for s in 1:length(selection)
        setlowerbound(prob.vars[:fs][s], 1)
    end

    status = solve(prob.model)
    if status == :Infeasible
        print_iis_gurobi(prob.model)
        error("Indicated infeasible on problem")
    end
    fullObjVal = getobjectivevalue(prob.model)

    # Setup warm start value
    for i in 1:length(prob.model.colVal)
        setvalue(Variable(prob.model, i), prob.model.colVal[i])
    end

    # Restoration phase since the bound change is in-place
    for b in 1:B
        for t in 1:T
            setlowerbound(prob.vars[:pg][b,t], orig_l_expand[b,t])
            setlowerbound(prob.vars[:h][b,t], orig_l_harden[b,t])
        end
    end

    for s in 1:length(selection)
        setlowerbound(prob.vars[:fs][s], 0)
    end

    return
end
