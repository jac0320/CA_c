function reactor(param::Dict, stoc::stocType, driver::Dict)

    pickIdx = rand([1:stoc.S;])
    pickScen = stoc.scenarios[pickIdx]

    expand = static_expansion(param, driver)
    closePlan = slr_shutdown(param, pickScen)
    harden = threshold_harden(param, pickScen, 0.0, true, closePlan)

    hDesign = designType(1, expand, harden)
    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info("[REACTOR] Cost is $totalcost = $expandcost + $hardencost")

    write_output_files(hDesign, driver)
    print_design(hDesign, param)

    return hDesign
end

function highland(param::Dict, stoc::stocType, driver::Dict)

    pickIdx = rand([1:stoc.S;])
    pickScen = stoc.scenarios[pickIdx]
    percSS = get_percentile_SS(stoc, 90, [pickIdx])

    expand = static_expansion(param, driver, hl=percSS)
    closePlan = slr_shutdown(param, pickScen)
    harden = threshold_harden(param, pickScen, 0.0, true, closePlan)

    hDesign = designType(1, expand, harden)
    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info("[HIGHLAND] The total cost is $totalcost = $expandcost + $hardencost")

    write_output_files(hDesign, driver)
    print_design(hDesign, param)

    return hDesign
end

function bathtub(param::Dict, stoc::stocType, driver::Dict)

    pickIdx = rand([1:stoc.S;])
    pickScen = stoc.scenarios[pickIdx]
    maxSS = get_percentile_SS(stoc, 100, [pickIdx])

    expand = static_expansion(param, driver)
    closePlan = slr_shutdown(param, pickScen)
    harden = threshold_harden(param, pickScen, maxSS, false, closePlan)

    hDesign = designType(1, expand, harden)
    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info(string("[BATHTUB]The total cost is $totalcost = $expandcost + $hardencost"))

    write_output_files(hDesign, driver)
    print_design(hDesign, param)

    return hDesign
end

function extreme(param::Dict, stoc::stocType, driver::Dict)

    pickIdx = rand([1:stoc.S;])
    pickScen = stoc.scenarios[pickIdx]
    maxSS = get_percentile_SS(stoc, 100, [pickIdx])

    expand = static_expansion(param, driver)
    closePlan = slr_shutdown(param, pickScen)
    harden = threshold_harden(param, pickScen, 1.1*maximum(pickScen.data["SS"]), false, closePlan)

    hDesign = designType(1, expand, harden)

    totalcost, expandcost, hardencost = get_design_cost(hDesign, param)
    info(string("[EXTREME]The total cost is $totalcost = $expandcost + $hardencost"))

    write_output_files(hDesign, driver)
    print_design(hDesign, param)

    return hDesign
end

# TODO this can be improved
function warmstart_heuristic(prob::oneProblem, param::Dict, stoc::stocType, driver::Dict; selection=[])

    isempty(selection) ? selection = [1:param[:S];] : selection = selection

    S, T, B = length(selection), param[:T], param[:B]

    expand = convert(Array{Int}, static_expansion(param, driver))

    orig_l_expand = convert(Array{Int}, getlowerbound(prob.vars[:pg]))
    orig_u_harden = convert(Array{Int}, getupperbound(prob.vars[:h]))
    orig_l_harden = convert(Array{Int}, getlowerbound(prob.vars[:h]))

    enforce_bound(prob, :pg, lb=expand)
    enforce_bound(prob, :h,  lb=orig_u_harden)
    enforce_bound(prob, :fs, lb=ones(Int, S))

    config_solver(prob.model, driver)
    status = solve(prob.model)
    status == :Infeasible && print_iis_gurobi(prob.model)
    status == :Infeasible && error("warm start failed")

    enforce_bound(prob, :pg, lb=orig_u_harden)
    enforce_bound(prob, :h,  lb=orig_l_expand)
    enforce_bound(prob, :fs, lb=zeros(Int, S))

    return
end
