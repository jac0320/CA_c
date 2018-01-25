##################################################################
#
#       Welcome to the Climate Adaptation project.
#       This is where everything starts.
#
##################################################################

module ClimateAdaptation

using JSON
using JuMP
using Iterators
using PowerModels
using Gurobi, MathProgBase  # Used for infeasible diagnostic

# myid() == 1 && config.PARALLEL && include("parallel.jl")

pkgdir = Pkg.dir("ClimateAdaptation")

include("$(pkgdir)/src/core/types.jl")
include("$(pkgdir)/src/core/input.jl")
include("$(pkgdir)/src/core/solver.jl")
include("$(pkgdir)/src/core/utility.jl")
include("$(pkgdir)/src/core/stoch.jl")
include("$(pkgdir)/src/core/soln.jl")

include("$(pkgdir)/src/formulation/general.jl")
# include("$(pkgdir)/src/formulation/evaluation.jl")
include("$(pkgdir)/src/formulation/sbd.jl")
include("$(pkgdir)/src/formulation/cuts.jl")

include("$(pkgdir)/src/algo/regular.jl")
include("$(pkgdir)/src/algo/shcg.jl")
include("$(pkgdir)/src/algo/shcgutility.jl")

# include("$(pkgdir)/src/algo/evaluation.jl")
# include("$(pkgdir)/src/algo/report.jl")
include("$(pkgdir)/src/algo/heuristic.jl")

export adcc

end # module
