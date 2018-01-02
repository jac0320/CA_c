function trim_bounds(prob, stoc, exargs; selection=[])

    isempty(selection) ? selection=[1:stoc.S;] : selection=selection

    protected_bus = []

    B = prob.param[:B]
    T = prob.param[:T]
    g_lim = zeros(B)
    h_u_lim = zeros(B)
    worst_ss = zeros(B)

    # Collect worst surge level at each time period
    for s in selection
        for t in 1:T
            for b in 1:B
                (stoc.scenarios[s].data["SS"][b,t] > worst_ss[b]) && (worst_ss[b] = stoc.scenarios[s].data["SS"][b,t])
            end
        end
    end

    # Using worst_ss conclude the upper bounds for h
    for b in 1:B
        h_u_lim[b] = ceil(max.(0.0, worst_ss[b] - prob.param[:Ele][b])/prob.param[:ProM][b])
    end

    # Trim the hardening upper bound with worst_ss
    for b in 1:B
        # if h_u_lim[b] > 0.0
        #     println("Regulating BUS $b hardeing upper bounds with $(h_u_lim[b])")
        # end
        for t in 1:T
            setupperbound(prob.vars[:h][b,t], h_u_lim[b])
        end
    end

    return prob
end
