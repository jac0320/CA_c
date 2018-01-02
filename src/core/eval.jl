"""
    This is a basic simulator for testing a design.
"""
function eval_basic(design::designType, simulator, formulation, prob::Dict, exargs::Dict;kwargs...)

    extras = Dict(kwargs)

    # Use simulator to generate scenarios
    stoc = simulator(B, S, T, exargs[stocMode])

    # Generate formulation
    results = fromulation(prob, param, stoc, design, exargs, "-")

    # Output results
    println("Evaluation results are ", results)

    return 0
end
