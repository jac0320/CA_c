function adcc(;kwargs...)

	driver = build_driver(Dict(kwargs))
	# *******************

    power = read_power(driver)
    stoc = get_scenarios(driver)
    param = read_parameters(power, stoc, driver)

	summary_driver_arguments(param, stoc, driver)
    summary_scenarios(stoc, param, driver)

    if driver[:ALGO] == "solve" || driver[:ALGO] == "regular"
        info(string("Sending original problem to solver ", driver[:ALGO]))
        totalTime = @elapsed problem, solution = deterministic(power, param, stoc, driver)
		println("Wall time [$totalTime]s")

    # elseif driver[:ALGO] == "sbd_heuristic" || driver[:ALGO] == "sbd"
    #     info(string("Running algorihtm with CG heuristic", driver[:ALGO]))
    #     totalTime = @elapsed problem, solution = sbd_heuristic(power, param, stoc, driver, sbd_master_formulation, sbd_subprob_formulation)
	# 	println("Wall time [$totalTime]s")
    #
    # elseif driver[:ALGO] == "sbd_norisk" || driver[:ALGO] == "sbdnr"
    #     info(string("Running algorithm : Sample-based Heuristic Decomposition SBD-NORISK. "))
    #     totalTime = @elapsed problem, solution = sbd_norisk(power, param, stoc, driver, sbd_subprob_formulation)
	# 	println("Wall time [$totalTime]s")
    #
	# elseif driver[:ALGO] == "heuristic" || driver[:ALGO] == "heu"
	# 	info(string("Running algorithm : Heuristic method ($(driver[:HEURISTIC])). "))
	# 	totalTime = @elapsed problem, solution = eval(parse(driver[:HEURISTIC]))(power, param, stoc, driver)
	# 	println("Wall time [$totalTime]s")
    #
	# elseif driver[:ALGO] == "damage_report"
	# 	info(string("Running reports : Reporting Stochastic Scenario Damages."))
	# 	analysis_scnearios(stoc, param)
    #
	# elseif driver[:ALGO] in ["solution_report", "sol_report", "sr"]
	# 	info("Running reports : solution summary.")
	# 	analysis_solution(power, param, stoc, driver)
    #
	# elseif driver[:ALGO] == "evaluate"
    #     info(string("Running algorihtms : Solution Evaluation."))
    #     evaluation(power, param, stoc, driver)
    #
    # elseif driver[:ALGO] == "enumerate"
	# 	info(string("Running enumertae algorithmic to explore the solution space"))
	# 	enumerator(power, param, stoc, driver)

	else
		error("Unkown algorithm (:ALGO) mode.")
    end

    info("\n\\\\------------ Successfully completed -------------//")
	return
end

function build_driver(args::Dict)

    driver = Dict()

	str2model = Dict("network"=>cnf_model,
					 "capacity"=>cb_model,
					 "dc"=>dcpf_model)

    haskey(args, :PROBLEM)      ? driver[:PROBLEM] = args[:T] : driver[:PROBLEM] = "ieee118"
    haskey(args, :T)            ? driver[:T] = args[:T] : driver[:T] = 5
    haskey(args, :EPS)          ? driver[:eps] = args[:EPS] : driver[:eps] = 0.0
    haskey(args, :STOCHFILE)    ? driver[:STOCHFILE] = args[:STOCHFILE] : driver[:STOCHFILE] = "test1.json"
	haskey(args, :STOCHMODE)	? driver[:STOCHMODE] = args[:STOCHMODE] : driver[:STOCHMODE] = ""
    haskey(args, :ALGO)         ? driver[:ALGO] = args[:ALGO] : driver[:ALGO] = "solve"
    haskey(args, :PARALLEL)     ? driver[:PARALLEL] = args[:PARALLEL] : driver[:PARALLEL] = false
    haskey(args, :CGHEURISTIC)  ? driver[:CGHEURISTIC] = args[:CGHEURISTIC] : driver[:CGHEURISTIC] = "improver_heu"
    haskey(args, :HEURISTIC)    ? driver[:HEURISTIC] = args[:HEURISTIC] : driver[:HEURISTIC] = "reactor"
    haskey(args, :EVALDESIGN)       ? driver[:EVALDESIGN] = args[:EVALDESIGN] : driver[:EVALDESIGN] = ""
    haskey(args, :EVALTARGET)       ? driver[:EVALTARGET] = args[:EVALTARGET] : driver[:EVALTARGET] = "feasibility"
    haskey(args, :COSTLambda)       ? driver[:COSTLambda] = args["COSTLambda"] : driver[:COSTLambda] = -1.0
    haskey(args, :SHEDLambda)       ? driver[:SHEDLambda] = args["SHEDLambda"] : driver[:SHEDLambda] = 0.95
    haskey(args, :CONGESTLambda)    ? driver[:CONGESTLambda] = args["CONGESTLambda"] : driver[:CONGESTLambda] = 1.0
    haskey(args, :ANGLESHIFTLambda) ? driver[:ANGLESHIFTLambda] = args["ANGLESHIFTLambda"] : driver[:ANGLESHIFTLambda] = 30.0
    haskey(args, :DISCOUNTLambda)   ? driver[:DISCOUNTLambda] = args["DISCOUNTLambda"] : driver[:DISCOUNTLambda] = -0.01
    haskey(args, :DEMANDLambda)     ? driver[:DEMANDLambda] = args["DEMANDLambda"] : driver[:DEMANDLambda] = 0.02
    haskey(args, :NAME)         ? driver[:NAME] = args[:NAME] : driver[:NAME] = "00"
    haskey(args, :WARMSTART)    ? driver[:WARMSTART] = args[:WARMSTART] : driver[:WARMSTART] = false

	haskey(args, :VERBOSE) 		? driver[:VERBOSE] = args[:VERBOSE] : driver[:VERBOSE] = false

	haskey(args, :MODEL) ? driver[:MODEL] = str2model[args[:MODEL]] : driver[:MODEL] = str2model["network"]
    haskey(args, :SOLVER)   ?  driver[:SOLVER] = args[:SOLVER] : driver[:SOLVER] = nothing

    haskey(args, :INPUTPATH) ? driver[:INPUTPATH] = args[:INPUTPATH] : driver[:INPUTPATH] = joinpath(homedir(),"Inputs","Climate")
    haskey(args, :OUTPUTPATH) ? driver[:OUTPUTPATH] = args[:OUTPUTPATH] : driver[:OUTPUTPATH] = joinpath(homedir(), "Outputs", "Climate")
	haskey(args, :TIMELIMIT) ? driver[:TIMELIMIT] = args[:TIMELIMIT] : driver[:TIMELIMIT] = 10800
    haskey(args, :TIMELIMITII) ? driver[:TIMELIMITII] = args[:TIMELIMITII] : driver[:TIMELIMITII] = 1800
	haskey(args, :TIMELIMITIII) ? driver[:TIMELIMITIII] = args[:TIMELIMITIII] : driver[:TIMELIMITIII] = 100
	haskey(args, :MAXITER) ? driver[:MAXITER] = args[:MAXITER] : driver[:MAXITER] = 999
	haskey(args, :SHOWLOG) ? driver[:SHOWLOG] = args[:SHOWLOG] : driver[:SHOWLOG] = false
	haskey(args, :BIGMINT) ? dirver[:BIGMINT] = args[:BIGMINT] : driver[:BIGMINT] = 9999999
	haskey(args, :BIGMDOUBLE) ? driver[:BIGMDOUBLE] = args[:BIGMDOUBLE] : driver[:BIGMDOUBLE] = 99999999.9999
	haskey(args, :WORKER) ? driver[:WORKER] = args[:WORKER] : driver[:WORKERS] = 2
	haskey(args, :OPTGAP) ? driver[:OPTGAP] = args[:OPTGAP] : driver[:OPTGAP] = 0.01
    haskey(args, :THREADS) ? driver[:THREADS] = args[:THREADS] : driver[:THREADS] = 8
    haskey(args, :WORKERTHREADS) ? driver[:WORKERTHREADS] = args[:WORKERTHREADS] : driver[:WORKERTHREADS] = 2
    haskey(args, :USESBDNORISK) ? driver[:USESBDNORISK] = args[:USESBDNORISK] : driver[:USESBDNORISK] = 0
    haskey(args, :WARMSTART) ? driver[:WARMSTART] = args[:USESBDNORISK] : driver[:WARMSTART] = true

    return driver
end
