"""
    Give an overview of the scenarios.
"""
function analysis_scnearios(stoc::stocType, param::Dict)
    S = stoc.S
    T = param[:T]
    B = param[:B]

    # Damage Report
    println(" ====== Overall surged bus (without any protection) ====== ")
    osb = Array{Int64}(S, T)
    fpc = Array{Any}(S, T)
    for s in 1:S
        timeline = "@ -> "
        costline = "@ -> "
        for t in 1:T
            sb = 0
            fpc_st = 0.0
            for b in 1:B
                if stoc.scenarios[s].data["SS"][b,t] > param[:Ele][b]
                    sb += 1
                    fpc_st += (stoc.scenarios[s].data["SS"][b,t] - param[:Ele][b]) / param[:ProM][b] * param[:Ch][b,t]
                end
             end
             osb[s,t] = sb
             fpc[s,t] = fpc_st
             timeline = string(timeline, osb[s,t], " -> ")
             costline = string(costline, round.(fpc[s,t],2), " -> ")
        end
        println("\n Scenarios $s : ", timeline, " !! $(maximum(osb[s, :]))")
		println("\t costline($s) : ", costline, " !! $(maximum(fpc[s, :]))")
    end
    println("=========================================================")
end

"""
    Give an overview of a solution
"""
function analysis_solution(power::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

    options = Dict(kwargs)

    S = stoc.S
    T = param[:T]
    B = param[:B]

    # outfname = "$(config.OUTPUTPATH)$(exargs[:OUTPATH])report_$(exargs[:NAME]).out"
    # outf = open(outfname, "w")

    if isfile(exargs[:EVALDESIGN])
        d = parse_design(exargs[:EVALDESIGN], param)
    elseif isfile("$(config.OUTPUTPATH)$(exargs[:EVALDESIGN])")
        d = parse_design("$(config.OUTPUTPATH)$(exargs[:EVALDESIGN])", param)
    else
        error("Cannot locate .json file")
    end
    pf = "REPORT: "

    print_design(d, param)

    # Summary Costs
    c_total, c_exp, c_har = get_design_cost(d, param)
    println("TOTAL-COST : $(c_total)")
    println("    EXPAND-COST : $(c_exp)")
    println("    HARDEN-COST : $(c_har)")

    println("--------------------")

    bt_exp = sum(d.pg[:,T])-sum(d.pg[:,1])
    bt_har = sum(d.h[:,T])-sum(d.h[:,1])
    bt_total = bt_exp + bt_har
    println("TOTAL-BUILD : $(bt_total)")
    println("    EXPAND-BUILD : $(bt_exp)")
    println("    HARDEN-BUILD : $(bt_har)")

    println("--------------------")

    # Calculate Build Average Height
    b_exp = d.pg[:,T] - param[:Pg0]
    b_har = d.h[:,T] - param[:H0]
    b_total = b_exp + b_har

    bh = Vector{Float64}()
    bh_exp = Vector{Float64}()
    bh_har = Vector{Float64}()
    for i in 1:B
        if b_total[i] > 0.0
            bh = [bh, repmat([param[:Ele][i]], b_total[i]);]
        end
        if b_exp[i] > 0.0
            bh_exp = [bh_exp, repmat([param[:Ele][i]], b_exp[i]);]
        end
        if b_har[i] > 0.0
            bh_har = [bh_har, repmat([param[:Ele][i]], b_har[i]);]
        end
    end

    baveh = mean(bh)
    baveh_exp = mean(bh_exp)
    baveh_har = mean(bh_har)

    bmedh = median(bh)
    bmedh_exp = median(bh_exp)
    bmedh_har = median(bh_har)

    bminh = minimum(bh)
    bminh_exp = minimum(bh_exp)
    bminh_har = minimum(bh_har)

    bmaxh = maximum(bh)
    bmaxh_exp = maximum(bh_exp)
    bmaxh_har = maximum(bh_har)

    bh_hist = fit(Histogram, bh, nbins=20, closed=:left)
    bh_exp_hist = fit(Histogram, bh_exp, 0:2.0:60.0, closed=:left)
    bh_har_hist = fit(Histogram, bh_har, 1.0:0.25:7.0, closed=:left)

    println("Average | Median | Max | Min-HEIGHT : $(baveh) | $(bmedh) | $(bmaxh) | $(bminh)")
    println("    EXPAND-Average | Median | Max | Min-HEIGHT : $(baveh_exp) | $(bmedh_exp) | $(bmaxh_exp) | $(bminh_exp)")
    println("    EXPAND-HIST : $(cumsum(reverse(bh_exp_hist.weights)))")
    println("    HARDEN-Average | Median | Max | Min-HEIGHT : $(baveh_har) | $(bmedh_har) | $(bmaxh_har) | $(bminh_har)")
    println("    HARDEN-HIST : $(cumsum(reverse(bh_har_hist.weights)))")

    println("--------------------")

    # Calculate Unit Costs Summary
    unit_exp = [sum((d.pg[:,i] - param[:Pg0])) for i in 1:T]
    unit_har = [sum((d.h[:,i] - param[:H0])) for i in 1:T]
    unit_mw_exp = [(d.pg[:,i] - param[:Pg0])'*param[:PgUB] for i in 1:T]
    unit_me_har = [(d.h[:,i] - param[:H0])'*param[:ProM] for i in 1:T]
    unit_cost_exp = [(d.pg[:,i] - param[:Pg0])'*param[:Cg][:,i] for i in 1:T]
    unit_cost_har = [(d.h[:,i] - param[:H0])'*param[:Ch][:,i] for i in 1:T]
    println("AVE-UNIT-EXPAND $(c_exp/unit_mw_exp[end])")
    println("AVE-UNIT-HARDEN $(c_har/unit_me_har[end])")
    println("UNIT-EXPAND $(unit_exp)")
    println("UNIT-HARDEN $(unit_har)")
    # println("BETA MW-EXPAND $(unit_mw_exp)")
    # println("BETA METER-HARDEN $(unit_me_har)")
    # println("BETA COST-EXPAND $(unit_cost_exp)")
    # println("BETA COST-HARDEN $(unit_cost_har)")
    # close(outf)
    return
end
