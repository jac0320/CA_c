function eval_basic(design::designType, simulator, formulation, prob::Dict, driver::Dict;kwargs...)

    extras = Dict(kwargs)
    stoc = simulator(B, S, T, driver[stocMode])
    results = fromulation(prob, param, stoc, design, driver, "-")
    println("Evaluation results are ", results)

    return 0
end
