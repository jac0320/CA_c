##################################################################
#
#       Welcome to the Climate Adaptation project.
#       This is where everything starts.
#
##################################################################

module ClimateAdaptation

using JSON
using ArgParse
using JuMP
using PowerModels
using StatsBase
using PowerModels
using Glob

# myid() == 1 && config.PARALLEL && include("parallel.jl")

pkgdir = Pkg.dir("ClimateAdaptation")

include("$(pkgdir)/src/types.jl")

include("$(pkgdir)/src/main.jl")
include("$(pkgdir)/src/core/config.jl")
include("$(pkgdir)/src/core/prob.jl")
include("$(pkgdir)/src/core/relax.jl")
include("$(pkgdir)/src/core/solver.jl")
include("$(pkgdir)/src/core/utility.jl")
include("$(pkgdir)/src/core/stoch.jl")
include("$(pkgdir)/src/core/soln.jl")

include("$(pkgdir)/src/formulation/general.jl")
include("$(pkgdir)/src/formulation/evaluation.jl")
include("$(pkgdir)/src/formulation/benders.jl")
include("$(pkgdir)/src/formulation/sbd.jl")
include("$(pkgdir)/src/formulation/cuts.jl")

include("$(pkgdir)/src/algo/deterministic.jl")
include("$(pkgdir)/src/algo/sbd_risk.jl")
include("$(pkgdir)/src/algo/sbd_norisk.jl")
include("$(pkgdir)/src/algo/sbd_utility.jl")
include("$(pkgdir)/src/algo/sbd_iso.jl")
include("$(pkgdir)/src/algo/sbd_heu.jl")
include("$(pkgdir)/src/algo/benders.jl")
include("$(pkgdir)/src/algo/evaluation.jl")
include("$(pkgdir)/src/algo/report.jl")
include("$(pkgdir)/src/algo/heuristic.jl")
include("$(pkgdir)/src/algo/enumerate.jl")

export adcc

end # module
