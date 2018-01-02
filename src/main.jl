using JSON
using ArgParse

include("core/config.jl")
config = read_config()

if myid() == 1 && config.PARALLEL
	include("parallel.jl")
end

using JuMP
using PowerModels
using Glob
using ProgressMeter
using StatsBase
using Compat
using DataFrames

include("core/types.jl")
using ADCCTypes	# Local Module

include("core/prob.jl")
include("core/relax.jl")
include("core/solver.jl")
include("core/utility.jl")
include("core/param.jl")
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

setup_envs()

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

    power = read_power(driver[:PROBLEM])

	# Driver is claimed to be the global driver
    T 	= driver[:T]
    S 	= driver[:S]
    driver[:B] = length(power["bus"])
	B	= driver[:B]

	# Generate or read stochasticity info into the problem
    stoc, driver = create_samples(B, S, T, driver, filePath = driver[:STOCHFILE])
    param, driver = get_parameters(power, stoc, driver)
	summary_driver_arguments(param, stoc, driver)
    summary_scenarios(stoc, param)

    if driver[:ALGO] == "det" || driver[:ALGO] == "regular"
        info(string("Sending original problem to solver ", driver[:ALGO]))
        totalTime = @elapsed problem, solution = deterministic(power, param, stoc, driver)
		println("Wall time [$totalTime]s")

    elseif driver[:ALGO] == "sbd_heuristic" || driver[:ALGO] == "sbd"
        info(string("Running algorihtm with CG heuristic", driver[:ALGO]))
        totalTime = @elapsed problem, solution = sbd_heuristic(power, param, stoc, driver, sbd_master_formulation, sbd_subprob_formulation)
		println("Wall time [$totalTime]s")

    elseif driver[:ALGO] == "sbd_norisk" || driver[:ALGO] == "sbdnr"
        info(string("Running algorithm : Sample-based Heuristic Decomposition SBD-NORISK. "))
        totalTime = @elapsed problem, solution = sbd_norisk(power, param, stoc, driver, sbd_subprob_formulation)
		println("Wall time [$totalTime]s")

	elseif driver[:ALGO] == "heuristic" || driver[:ALGO] == "heu"
		info(string("Running algorithm : Heuristic method ($(driver[:HEURISTIC])). "))
		totalTime = @elapsed problem, solution = eval(parse(driver[:HEURISTIC]))(power, param, stoc, driver)
		println("Wall time [$totalTime]s")

	elseif driver[:ALGO] == "damage_report"
		info(string("Running reports : Reporting Stochastic Scenario Damages."))
		analysis_scnearios(stoc, param)

	elseif driver[:ALGO] in ["solution_report", "sol_report", "sr"]
		info("Running reports : solution summary.")
		analysis_solution(power, param, stoc, driver)

    elseif driver[:ALGO] == "resample"
		info(string("Resampling from a large set of senarios"))
		for i in 1:driver[:REPEAT]
			srand(i)
			randpick = randperm(stoc.S)
			info("Picking scenario subset [$(driver[:SUBSETS])] $(randpick[1:driver[:SUBSETS]])")
			substoc = subsetting_stocType(stoc, randpick[1:driver[:SUBSETS]], driver)
			info("Writing stoc file into output folder... NAME=$(driver[:NAME])")
			write_stocType_json(substoc, string(driver[:NAME],"_$i.json"))
		end

	elseif driver[:ALGO] == "evaluate"
        info(string("Running algorihtms : Solution Evaluation."))
        evaluation(power, param, stoc, driver)

    elseif driver[:ALGO] == "enumerate"
		info(string("Running enumertae algorithmic to explore the solution space"))
		enumerator(power, param, stoc, driver)

	elseif driver[:ALGO] == "benders"		# Not finished
		error("This feature is not ready...")
        info(string("Running algorihtm ", driver[:ALGO]))
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
    end

    info("\n\\\\------------ Successfully completed -------------//")
	return
end

function summary_driver_arguments(param::Dict, stoc::stocType, driver::Dict)

    info("Problem Instance      : ", driver[:PROBLEM])
    info("Characteristic        : ", driver[:MODEL])
    info("Stochastic Mode       : ", driver[:STOCHMODE])
    info("Algorithm             : ", driver[:ALGO])
    info("Time Periods(T)       : ", driver[:T])
    info("Scenario Count(S)     : ", driver[:S])
    info("Risk(eps)             : ", driver[:eps])
    info("Features              : ", driver[:FEATURES])
	info("Demand Change         : ", driver[:DEMANDLambda])
    info("Shedding Allowing     : ", driver[:SHEDLambda])
	info("Congestion            : ", driver[:CONGESTLambda])
	info("Angle Shift Limit     : ", driver[:ANGLESHIFTLambda])
	info("Discounting Cost      : ", driver[:DISCOUNTLambda])
	info("Cost Ratio (expand)   : ", driver[:COSTLambda])
	info("Time Limit L1         : ", config.TIMELIMIT)
	info("Time Limit L2         : ", config.TIMELIMITII)
	info("Time Limit L3         : ", config.TIMELIMITIII)
	info("Parallel Indicator    : ", config.PARALLEL)
	info("Workers Utilized      : ", config.WORKERS)
    info("Single Worker Threads : ", config.WORKERTHREADS)
	info("Workers Count on Node : ", config.JOBPERWORKER)
	info("Warm Start            : ", config.WARMSTART)
	info("Job created at        : ", now())
	info("Job Output Name       : ", driver[:NAME])
	info("Julia Version         : ", VERSION)
	info("CPU Cores             : ", Sys.CPU_CORES)
	info("Machine INFO          : ", Sys.MACHINE)
	info("CPU Summary           : ")
	Sys.cpu_summary()
	info("Head Node Name        : ")
	run(`hostname`)
	info("ADDITIONAL CSV		: ", !(driver[:CSV] == nothing))

    @assert param[:S]   == stoc.S
    @assert param[:S]   == driver[:S]
    @assert param[:T]   == stoc.T
    @assert param[:T]   == driver[:T]
    @assert param[:eps] == driver[:eps]

	return
end
