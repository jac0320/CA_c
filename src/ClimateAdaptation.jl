##################################################################
#
#       Welcome to the Climate Adaptation project.
#       This is where everything starts.
#
##################################################################

module ClimateAdaptation

using JSON
using ArgParse

include("core/config.jl")
config = read_config()

# myid() == 1 && config.PARALLEL && include("parallel.jl")

using JuMP, PowerModels, StatsBase
using PowerModels
using Glob, ProgressMeter
using DataFrames

# Enterance, can use : commandline arguments | functional arguments | Default setting
include("$(Pkg.dir("ClimateAdaptation"))/src/main.jl")

export adcc

end # module
