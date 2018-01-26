using JSON
using ArgParse

include("core/config.jl")
config = read_config()

if myid() == 1 && config.PARALLEL
	include("parallel.jl")
end

using JuMP, PowerModels, StatsBase
using PowerModels
using Glob, ProgressMeter
using DataFrames

include("core/types.jl")
using ADCCTypes	# Local Module

include("core/prob.jl")
include("core/relax.jl")
include("core/solver.jl")
include("core/utility.jl")
include("core/stoch.jl")
include("core/soln.jl")

include("formulation/general.jl")
include("formulation/evaluation.jl")
include("formulation/benders.jl")
include("formulation/sbd.jl")
include("formulation/cuts.jl")

include("algo/deterministic.jl")
include("algo/sbd_risk.jl")
include("algo/sbd_norisk.jl")
include("algo/sbd_utility.jl")
include("algo/sbd_iso.jl")
include("algo/sbd_heu.jl")
include("algo/benders.jl")
include("algo/evaluation.jl")
include("algo/report.jl")
include("algo/heuristic.jl")
include("algo/enumerate.jl")


if config.SOLVER == "Gurobi" || config.SOLVER == "gurobi" || config.SOLVER=="GUROBI"
	using Gurobi
elseif config.SOLVER == "CPLEX" || config.SOLVER == "Cplex" || config.SOLVER == "cplex"
	using CPLEX
elseif config.SOLVER == "Cbc" || config.SOLVER == "cbc"
	using Cbc
end

function adcc(;kwargs...)

	# ======================================================================== #
	inputArgs = Dict(kwargs)
	args = parse_commandline_args(inputArgs)
	driver = get_driver_args(
	    args,	# Command line input
	    PROBLEM = "nesta_case14_ieee.m",
	    MODEL = network_characteristic,
	    STOCHMODE = "evolving",
	    STOCHFILE = "./",
	    ALG = "regular",
	    evalDesign = "./",
	    evalTarget = "feasibility",
		FEATURES = ["sample-based-risk","surge-load-shed"],	# FINAL DECISION
	 	EPS= 0.0,
	    T = 6,
	    S = 10
	)
	# ======================================================================== #

    power = read_power(driver)
    stoc = get_scenarios(driver)
    param = read_parameters(power, stoc, driver)

	summary_driver_arguments(param, stoc, driver)
    summary_scenarios(stoc, param)

    if driver[:ALGO] == "solve" || driver[:ALGO] == "regular"
        println(string("Sending original problem to solver ", driver[:ALGO]))
        totalTime = @elapsed problem, solution = deterministic(power, param, stoc, driver)
		println("Wall time [$totalTime]s")

    elseif driver[:ALGO] == "sbd_heuristic" || driver[:ALGO] == "sbd"
        println(string("Running algorihtm with CG heuristic", driver[:ALGO]))
        totalTime = @elapsed problem, solution = sbd_heuristic(power, param, stoc, driver, sbd_master_formulation, sbd_subprob_formulation)
		println("Wall time [$totalTime]s")

    elseif driver[:ALGO] == "sbd_norisk" || driver[:ALGO] == "sbdnr"
        println(string("Running algorithm : Sample-based Heuristic Decomposition SBD-NORISK. "))
        totalTime = @elapsed problem, solution = sbd_norisk(power, param, stoc, driver, sbd_subprob_formulation)
		println("Wall time [$totalTime]s")

	elseif driver[:ALGO] == "heuristic" || driver[:ALGO] == "heu"
		println(string("Running algorithm : Heuristic method ($(driver[:HEURISTIC])). "))
		totalTime = @elapsed problem, solution = eval(parse(driver[:HEURISTIC]))(power, param, stoc, driver)
		println("Wall time [$totalTime]s")

	elseif driver[:ALGO] == "damage_report"
		println(string("Running reports : Reporting Stochastic Scenario Damages."))
		analysis_scnearios(stoc, param)

	elseif driver[:ALGO] in ["solution_report", "sol_report", "sr"]
		println("Running reports : solution summary.")
		analysis_solution(power, param, stoc, driver)

    elseif driver[:ALGO] == "resample"
		println(string("Resampling from a large set of senarios"))
		for i in 1:driver[:REPEAT]
			srand(i)
			randpick = randperm(stoc.S)
			println("Picking scenario subset [$(driver[:SUBSETS])] $(randpick[1:driver[:SUBSETS]])")
			substoc = subsetting_stocType(stoc, randpick[1:driver[:SUBSETS]], driver)
			println("Writing stoc file into output folder... NAME=$(driver[:NAME])")
			write_stocType_json(substoc, string(driver[:NAME],"_$i.json"))
		end

	elseif driver[:ALGO] == "evaluate"
        println(string("Running algorihtms : Solution Evaluation."))
        evaluation(power, param, stoc, driver)

    elseif driver[:ALGO] == "enumerate"
		println(string("Running enumertae algorithmic to explore the solution space"))
		enumerator(power, param, stoc, driver)

	elseif driver[:ALGO] == "benders"		# Not finished
		error("This feature is not ready...")
        println(string("Running algorihtm ", driver[:ALGO]))
		benders(climate_benders_master_variables,
				climate_benders_master_constraints,
				climate_bender_master_objective,
				climate_benders_subprob_variables,
				climate_benders_subprob_constraints,
				climate_benders_subprob_objective,
				nothing, 	# Characteristic formulation is currently not supported
				power=power,
				param=param,
				stoc=stoc,
				exargs=driver)

	elseif driver[:ALGO] == "simulation"
		error("No implementation yet about this algorithm.")

	elseif driver[:ALGO] == "generate"
		error("No implementation yet about generting random scenarios.")
    end

    println("\n\\\\------------ Successfully completed -------------//")
	return
end

function summary_driver_arguments(param::Dict, stoc::stocType, driver::Dict)

    println("Problem Instance      : ", driver[:PROBLEM])
    println("Characteristic        : ", driver[:MODEL])
    println("Stochastic Mode       : ", driver[:STOCHMODE])
    println("Algorithm             : ", driver[:ALGO])
    println("Time Periods(T)       : ", driver[:T])
    println("Scenario Count(S)     : ", driver[:S])
    println("Risk(eps)             : ", driver[:eps])
    println("Features              : ", driver[:FEATURES])
	println("Demand Change         : ", driver[:DEMANDLambda])
    println("Shedding Allowing     : ", driver[:SHEDLambda])
	println("Congestion            : ", driver[:CONGESTLambda])
	println("Angle Shift Limit     : ", driver[:ANGLESHIFTLambda])
	println("Discounting Cost      : ", driver[:DISCOUNTLambda])
	println("Cost Ratio (expand)   : ", driver[:COSTLambda])
	println("Time Limit L1         : ", config.TIMELIMIT)
	println("Time Limit L2         : ", config.TIMELIMITII)
	println("Time Limit L3         : ", config.TIMELIMITIII)
	println("Parallel Indicator    : ", config.PARALLEL)
	println("Workers Utilized      : ", config.WORKERS)
    println("Single Worker Threads : ", config.WORKERTHREADS)
	println("Workers Count on Node : ", config.JOBPERWORKER)
	println("Warm Start            : ", config.WARMSTART)
	println("Job created at        : ", now())
	println("Job Output Name       : ", driver[:NAME])
	println("Julia Version         : ", VERSION)
	println("CPU Cores             : ", Sys.CPU_CORES)
	println("Machine INFO          : ", Sys.MACHINE)
	println("CPU Summary           : ")
	Sys.cpu_summary()
	println("Head Node Name        : ")
	run(`hostname`)

    @assert param[:S]   == stoc.S
	@assert param[:B]   == driver[:B]
    @assert param[:S]   == driver[:S]
    @assert param[:T]   == stoc.T
    @assert param[:T]   == driver[:T]
    @assert param[:eps] == driver[:eps]

	return
end
