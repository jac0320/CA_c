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
    #     totalTime = @elapsed problem, solution = sbd_heuristic(power, param, stoc, driver, sbd_master_formulation, build_sp)
	# 	println("Wall time [$totalTime]s")
    #
    # elseif driver[:ALGO] == "sbd_norisk" || driver[:ALGO] == "sbdnr"
    #     info(string("Running algorithm : Sample-based Heuristic Decomposition SBD-NORISK. "))
    #     totalTime = @elapsed problem, solution = sbd_norisk(power, param, stoc, driver, build_sp)
	# 	println("Wall time [$totalTime]s")
    #
	
	elseif driver[:ALGO] == "heuristic" || driver[:ALGO] == "heu"
		info(string("Running algorithm : Heuristic method ($(driver[:HEURISTIC])). "))
		totalTime = @elapsed solution = eval(parse(driver[:HEURISTIC]))(power, param, stoc, driver)
		println("Wall time [$totalTime]s")
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

    elseif driver[:ALGO] == "enumerate"
		info(string("Running enumertae algorithmic to explore the solution space"))
		enumerator(power, param, stoc, driver)
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

function read_power(exargs::Dict)

	pname = exargs[:PROBLEM]
	if pname in ["ieee14", "ieee118", "ieee118p", "gs4", "lmbd3"]
		power = PowerModels.parse_file(joinpath(Pkg.dir("ClimateAdaptation"), "instance", "$(pname).m"))
		info("Finish reading $(pname) power data...")
		return power
	else
		error("No instance network found.")
	end
end

function read_parameters(power::Dict, stoc::stocType, exargs::Dict)

    if isfile(joinpath(Pkg.dir("ClimateAdaptation"),"instance","$(exargs[:PROBLEM])_param.json"))
        pF = JSON.parsefile(joinpath(Pkg.dir("ClimateAdaptation"),"instance","$(exargs[:PROBLEM])_param.json"))
    else
        error("ERROR|param.jl|import_param()|No parameter file detected.")
    end

    p = init_parameters(power, stoc, exargs)
    p = populate_basic(p, pF)
    p = populate_costs(p, exargs)
    p = populate_g0(p, exargs)
    p = populate_demand(p, exargs)
    p = populate_edge(p, exargs)
    p = get_aslDet(p, stoc)
    p = get_assDet(p, stoc)
    p = resolve_undersea_load_shift(p, stoc)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    p = check_parameter_intactness(p)

	sync_parameters(p, exargs)

    return p
end

function init_parameters(power::Dict, stoc::stocType, exargs::Dict)

    # Parameter Dictionary Bundle
    param = Dict()
    bus = Dict()
    gen = Dict()
    line = Dict()

    for bidx in keys(power["bus"])
        bus[Int(parse(bidx))] = power["bus"][bidx]
    end

    for gidx in keys(power["gen"])
        gen[Int(parse(gidx))] = power["gen"][gidx]
    end

    for lidx in keys(power["branch"])
        line[Int(parse(lidx))] = power["branch"][lidx]
    end

    line = filter((i, branch) -> branch["br_status"] == 1, line)

    # Dimension defining parameters
    param[:bus]     = bus
    param[:gen]     = gen
    param[:line]    = line
    param[:S]       = S   = stoc.S
    param[:T]       = T   = exargs[:T]

    correctEps = exargs[:eps]
    exargs[:eps] = correctEps
    param[:eps]     = exargs[:eps]

    param[:B]       = B = bus.count;
    param[:L]       = L = line.count;
    param[:G]       = gen.count;
    param[:Cg]      = zeros(Float64, B, T)   # Costs changes over time
    param[:Ch]      = zeros(Float64, B, T)
    param[:Pd] 		= zeros(Float64, B, T)	# Power Demand :: ?Power data || uncertainty
    param[:PgUB]    = zeros(Float64, B, T)
    param[:Pg0]     = zeros(Int, B)
    param[:H0]      = zeros(Int, B)
    param[:XCoord]  = zeros(Float64, B)
    param[:YCoord]  = zeros(Float64, B)
    param[:Ele]     = zeros(Float64, B)
    param[:ProM]    = zeros(Float64, B)
    param[:Pgbar] 	= zeros(Float64, B)		# Upper bound of the building generator limits
    param[:Hbar]    = zeros(Float64, B)
    param[:EDGE]    = zeros(Int, B, B)		# 0/1 Matrix to indicate the line condition
    param[:Edge]    = Dict()
    param[:Lcap]	= zeros(Float64, B, B)	# Line Capacity :: ?Power datasol
	param[:lineX]	= zeros(Float64, B, B)
    param[:AngleLimit] = pi;
    param[:AngleShiftLimit] = pi / 180 * exargs[:ANGLESHIFTLambda];
    param[:SHEDLambda] = exargs[:SHEDLambda];
    param[:surplusLoad] = 0.0

    return param
end

function populate_basic(param::Dict, pF::Dict)

    B = param[:B]

    haskey(pF, "XCoord") ? param[:XCoord] = copy(pF["XCoord"]) : warn("Parameter XCoord not found")
    haskey(pF, "YCoord") ? param[:YCoord] = copy(pF["YCoord"]) : warn("Parameter YCoord not found")
    haskey(pF, "Ele") ? param[:Ele] = copy(pF["Ele"]) : error("Parameter Ele missing")
    haskey(pF, "PgUB") ? param[:PgUB] = copy(pF["PgUB"]) : error("Parameter PgUB missing")
    haskey(pF, "Cg") ? param[:CgOrig] = copy(pF["Cg"]) : error("Parameter Cg missing")
    haskey(pF, "Ch") ? param[:ChOrig] = copy(pF["Ch"]) : error("Parameter Ch missing")
    haskey(pF, "ProM") ? param[:ProM] = copy(pF["ProM"]) : error("Parameter ProM missing")
    haskey(pF, "Pgbar") ? param[:Pgbar] = copy(pF["Pgbar"]) : error("Parameter Pgbar missing")
    haskey(pF, "Hbar") ? param[:Hbar] = copy(pF["Hbar"]) : error("Parameter Hbar missing")
    if haskey(pF, "RefBus")
        param[:RefBus] = pF["RefBus"]
    else
        param[:RefBus] = rand(1:B)
        warn("Parameter RefBus not found, using random slackbus")
    end
    return param
end

function check_parameter_intactness(param::Dict)

    B = param[:B]
    T = param[:T]

	if length(param[:Cg]) == B
		temp = [];
		for i in 1:T temp=[temp;param[:Cg]] end
		param[:Cg] = temp;
	elseif length(param[:Cg]) == B * T
		param[:Cg] = param[:Cg];
	else
		error("ERROR|model.jl|param[:Ch]|check_parameter_intactness()::Parameters missing(param[:Cg]).");
	end
	param[:Cg] = reshape(param[:Cg], B, T);	# Time in-variant costs structure

	if length(param[:Ch]) == B
		temp = []
		for i in 1:T temp=[temp;param[:Ch]] end
		param[:Ch] = temp;
	elseif length(param[:Ch]) == B * T
		param[:Ch] = param[:Ch];
	else
		error("ERROR::model.jl::param[:Ch]eck_parameter_intactness()::Parameters missing(param[:Ch]).");
	end
	param[:Ch] = reshape(param[:Ch], B, T)

	return param
end

function populate_edge(param::Dict, exargs::Dict)

    B = param[:B]
    L = param[:L]

    for b in 1:B
        param[:Pd][b,1] = param[:bus][b]["pd"] * 100;
        param[:Edge][b] = Dict();
        param[:Edge][b]["out"] = [];
        param[:Edge][b]["in"] = [];
    end

    for l in 1:L
        from = param[:line][l]["f_bus"]
        to = param[:line][l]["t_bus"]
        param[:Lcap][from, to] = param[:line][l]["rate_a"] * 100 * exargs[:CONGESTLambda]
        param[:EDGE][from, to] = 1
        if (from in param[:Edge][to]["in"])
            # info("Found double coordior ($to, $from)")
            param[:Lcap][from, to] += param[:Lcap][from, to]
        else
            push!(param[:Edge][to]["in"],from)
            push!(param[:Edge][from]["out"],to)
        end
    end

	E = length(param[:line])
	for l in 1:E
		param[:lineX][param[:line][l]["f_bus"], param[:line][l]["t_bus"]] = param[:line][l]["br_x"] / 100
	end

    return param
end

function populate_demand(param::Dict, exargs::Dict)

    T = param[:T]
    B = param[:B]

    for b in 1:B
        param[:Pd][b,1] = param[:bus][b]["pd"] * 100;
    end
    for t in 2:T
        for b in 1:B
            param[:Pd][b, t] = param[:Pd][b, t-1]*(1+exargs[:DEMANDLambda])
        end
    end
    info("Total Demand Growth $(sum(param[:Pd][:,T])-sum(param[:Pd][:,1])) MW")

    return param
end

function populate_costs(param::Dict, exargs::Dict)

	Cg = param[:CgOrig]
	Ch = param[:ChOrig]

    Cg, Ch = lambdaCostTune(Cg, Ch, exargs[:COSTLambda])

    B = param[:B]
    T = param[:T]

    for b in 1:B
        for t in 1:T
            param[:Cg][b,t] = Cg[b] * (1 + exargs[:DISCOUNTLambda]^(t-1))
            param[:Ch][b,t] = Ch[b] * (1 + exargs[:DISCOUNTLambda]^(t-1))
        end
    end

    return param
end

function populate_g0(param::Dict, exargs::Dict)

    G = param[:G]
    for g in 1:param[:G]
        if param[:gen][g]["pg"] > 0.1
            param[:Pg0][param[:gen][g]["gen_bus"]] = 1.0
        end
    end
    return param
end

function get_aslDet(param::Dict, stoc::stocType)

    B = param[:B]
    T = param[:T]
    S = param[:S]

    aslDET = zeros(Bool, B, T, S)
    param[:Ele]
    for s in 1:S
        for i in 1:B
            for t in 1:T
                if stoc.scenarios[s].data["SL"][t] > param[:Ele][i]
                    # @show "SLR Flooding Scenario $(s) :: BUS $(i) || TIME $(t)  (ELE $(param[:Ele][i]), SLR $(stoc.scenarios[s].data["SL"][t]))"
                    aslDET[i,t,s] = false
                else
                    aslDET[i,t,s] = true
                end
            end
        end
    end

    param[:aslDet] = aslDET

    return param
end

function get_assDet(param::Dict, stoc::stocType)

    B = param[:B]
    T = param[:T]
    S = param[:S]

    assDET = zeros(Bool, B, T, S)

    for s in 1:S
        for i in 1:B
            for t in 1:T
                if stoc.scenarios[s].data["SS"][i,t] >= param[:Ele][i]
                    assDET[i,t,s] = false
                else
                    assDET[i,t,s] = true
                end
            end
        end
    end

    param[:assDet] = assDET

    return param
end

function lambdaCostTune(Cg::Array, Ch::Array, lambda::Float64)

    # Only deals with one time period (initial time period)
    if lambda < 0    # Can be disabled
        return Cg, Ch
    else

        B = length(Cg)

        totalCg = sum(Cg)
        totalCh = sum(Ch)
        origLambda = totalCg / (totalCg+totalCh)

        info("Original lambda $origLambda vs $lambda")

        totalC = totalCg + totalCh

        CgRatio = Array{Float64}(length(Cg))
        ChRatio = Array{Float64}(length(Ch))

        CgRatio = Cg/totalCg
        ChRatio = Ch/totalCh

        CgLambda = lambda * totalC * CgRatio
        ChLambda = (1-lambda) * totalC * ChRatio

        return CgLambda, ChLambda

    end
end

function resolve_undersea_load_shift(param::Dict, stoc::stocType)

    origLoad = copy(param[:Pd])

    B = param[:B]
    T = param[:T]
    S = stoc.S

    # Risk Parameter Validation
    @assert S == param[:S]

    feaLoadShift = true

    @assert haskey(param, :aslDet) # Will not proceed if this cannot go pass

    aslDet = param[:aslDet]

    param[:aslDetPd] = zeros(B, T, S)

    for i in 1:B
        SHIFT = zeros(Float64, T, S)
        for s in 1:S
            for t in 1:T
                nom = 0.0
                denom = 0.0
                for j in 1:B
                    # All loss load
                    nom += param[:Pd][j,t] * (1-aslDet[j, t, s])
                    # All operational node's load
                    denom += aslDet[j,t,s]*param[:Pd][j,t]
                end
                SHIFT[t,s] = nom/denom
            end
        end

        for t in 1:T
            for s in 1:S
                param[:aslDetPd][i,t,s] = origLoad[i,t] - origLoad[i,t]*(1-param[:aslDet][i,t,s]) + param[:aslDet][i,t,s]*origLoad[i,t]*SHIFT[t,s]
                if param[:aslDetPd][i,t,s] < 0
                    warn("Issue with bus $i, time point $t, and sample $s.")
                end
            end
        end
    end

    return param
end

function sync_parameters(param::Dict, exargs::Dict)
	exargs[:B] = param[:B]
end
