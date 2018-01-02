type configType
	INPUTPATH 	   ::AbstractString
	OUTPUTPATH	   ::AbstractString
	TIMELIMIT	   ::Float64
	TIMELIMITII	   ::Float64
	TIMELIMITIII   ::Float64
	TIMELIMITIV    ::Float64
	SOLVER		   ::AbstractString
    ENVS           ::Any
	RUNSEED		   ::Int
	TOLERANCE	   ::Float64
	MAXITER        ::Int
	SHOWLOG		   ::Int
	SCREENSHOW	   ::Int
	BIGMINT		   ::Int
	BIGMDOUBLE	   ::Float64
	FILEOUTPUT	   ::Int
	CONVERGENCE    ::Float64
	VERBOSE		   ::Int
	PARALLEL       ::Bool
	WORKERS		   ::Array
    JOBPERWORKER   ::Int
	OPTGAP		   ::Float64
	THREADS		   ::Int
    MAINTHREADS    ::Int
    WORKERTHREADS  ::Int
    USESBDNORISK   ::Bool
	WARMSTART	   ::Bool
	configType() = new()
end

function read_config(configPath::AbstractString="")

    config = configType()

	userARGS = parse_commandline_args()

    # Initialization of the configurations
    config.INPUTPATH = pwd()"/instance/"
    config.OUTPUTPATH = pwd()"/instance/"
    config.TIMELIMIT = 300
	config.TIMELIMITII = 100
	config.TIMELIMITIII = 50
	config.TIMELIMITIV = 10
    config.SOLVER = "Cbc"
    config.RUNSEED = 9812739127391
    config.TOLERANCE = 0.0001
	config.MAXITER = 999
    config.SCREENSHOW = 0
    config.BIGMINT = 9999999
    config.BIGMDOUBLE = 99999999.9999
    config.FILEOUTPUT = 0
    config.CONVERGENCE = 0.001
    config.VERBOSE = 0
    config.PARALLEL = 0
    config.WORKERS = [];   # Active when parallel structure is turned on.
	config.JOBPERWORKER = 1
    config.OPTGAP = 0.01
	config.THREADS = 0
	config.MAINTHREADS = 16
	config.WORKERTHREADS = 4
	config.USESBDNORISK = 0
	config.WARMSTART = 0

    # Read configuration from the config.txt file
	localRepo = false
	currentRepo = false
	if isfile(homedir()"/Debug/Climate/config.json")
		configPath = homedir()"/Debug/Climate/config.json"
		cfd = JSON.parsefile(configPath)
		info(string("Configuration directory : ", configPath))
		localRepo = true
	else
		info("No configuration file detected. Using default setting.")
		return config
	end

	if localRepo || currentRepo
	    config.INPUTPATH = cfd["INPUTPATH"]
	    config.OUTPUTPATH = cfd["OUTPUTPATH"]
	    config.TIMELIMIT = cfd["TIMELIMIT"]
		config.TIMELIMITII = cfd["TIMELIMITII"]
		config.TIMELIMITIII = cfd["TIMELIMITIII"]
		config.TIMELIMITIV = cfd["TIMELIMITIV"]
	    config.RUNSEED = cfd["RUNSEED"]
	    config.SOLVER = cfd["SOLVER"]
	    config.TOLERANCE = cfd["TOLERANCE"]
		config.MAXITER = cfd["MAXITER"]
	    config.SHOWLOG = cfd["SHOWLOG"]
	    config.SCREENSHOW = cfd["SCREENSHOW"]
	    config.BIGMINT = cfd["BIGMINT"]
	    config.BIGMDOUBLE = cfd["BIGMDOUBLE"]
	    config.FILEOUTPUT = cfd["FILEOUTPUT"]
	    config.CONVERGENCE = cfd["CONVERGENCE"]
	    config.VERBOSE = cfd["VERBOSE"]
	    config.PARALLEL = cfd["PARALLEL"]
	    config.WORKERS = cfd["WORKERS"]
		config.JOBPERWORKER = cfd["JOBPERWORKER"]
	    config.OPTGAP = cfd["OPTGAP"]
		config.THREADS = cfd["THREADS"]
		config.MAINTHREADS = cfd["MAINTHREADS"]
		config.WORKERTHREADS = cfd["WORKERTHREADS"]
		config.USESBDNORISK = cfd["USESBDNORISK"]
		config.WARMSTART = cfd["WARMSTART"]
	else
		info("No configuration file detected. Taking default value.")
	end


	(userARGS["SOLVER"] != "config") && (config.SOLVER = userARGS["SOLVER"])
	(userARGS["PARALLEL"] != "config") && (config.PARALLEL = Bool(parse(userARGS["PARALLEL"])))
	(userARGS["USESBDNORISK"] != "config") && (config.USESBDNORISK = Bool(parse(userARGS["USESBDNORISK"])))
	(userARGS["WORKERS"] > 0) && (config.WORKERS = [userARGS["WORKERS"]])
	(userARGS["WARMSTART"] != "config") && (config.WARMSTART = Bool(parse(userARGS["WARMSTART"])))
	(userARGS["TIMELIMIT"] > 0.0) && (config.TIMELIMIT = userARGS["TIMELIMIT"])
	(userARGS["TIMELIMITII"] > 0.0) && (config.TIMELIMITII = userARGS["TIMELIMITII"])
	(userARGS["MAXITER"] > 0) && (config.TIMELIMITIII = userARGS["MAXITER"])
	(userARGS["SEED"] != "config") && (config.RUNSEED = Int(parse(userARGS["SEED"])))

    return config
end

"""
    Parse commandline arguments input the code.
    If users prefer to call using functional arguments.
    It will do the same work.
"""
function parse_commandline_args(inputDriver::Dict=Dict())

    s = ArgParseSettings()

	@add_arg_table s begin
        "PROBLEM"
            help = "A problem/network identifier or path towards the .m file. (| ieee14 | nesta_case14_ieee |)"
            arg_type = AbstractString
            default = "ieee14"
        "MODEL"
            help = "Model type. (| capacity | network | dc | ac (not yet) |)"
            arg_type = AbstractString
            default = "network"
        "ALG"
            help = "Choice of algorithm to go with. (enumerate | regular | evaluate | sbd_heuristic | sbd_norisk | heuristic | benders (not yet) |)"
            arg_type = AbstractString
            default = "regular"
        "--T"
            help = "Total time steps"
            arg_type = Int
            default = 0
        "--S"
            help = "Target sample counts for multiple usage. (For internal generated scenarios only.)"
            arg_type = Int
            default = 0
        "--EPS"
            help = "Choice of risk allowed for the problem to be tested. If >= 1, then treat as counts for allowed infeasible scenarios."
            arg_type = Float64
            default = 0.0
        "--STOCHMODE", "--sm"
            help = "Internally generate samples with a specific mode. (| File | peaceful | evolving | random | fierce | evolving-highslr | evolving-fierce | evolving-95-perc | evolving-average | evolving-highslr-average | evolving-fierce-average |)"
            arg_type = AbstractString
            default = "evolving"
        "--STOCHFILE", "--sf"
            help = "User specified stoch file (.json) for reading in. This overrides the internal scenario generation commands."
            arg_type = AbstractString
            default = ""
        "--SOLVER"
            help = "User specified choice of solvers to use. (overrise configuration file)"
            arg_type = AbstractString
            default = "config"
        "--PARALLEL", "--PAR"
            help = "Enable parallel structure or not. (Override configuration file)"
            arg_type = AbstractString
            default = "config"
		"--WARMSTART", "--ws"
			help = "Apply warm start scheme at several specific locations. (Override configuration file)"
			arg_type = AbstractString
			default = "config"
		"--WORKERS"
			help = "Specify number workers would like to utilize for solving. (Override configuration file)"
			arg_type = Int
			default = -1
		"--TIMELIMIT", "--TL1"
			help = "Sepcify the total solver running time for level I control. (Override configuration file)"
			arg_type = Float64
			default = -1.0
		"--TIMELIMITII", "--TL2"
			help = "Specify the total solver running time for level II control. (Override configuration file)"
			arg_type = Float64
			default = -1.0
		"--MAXITER"
			help = "Specify the total iterations when certain heuristic algorithm is utilized. (Override configuration file)"
			arg_type = Int
			default = -1
		"--LEVEL"
			help = "Enumertae solution depth. (No more than 3 is recommanded)"
			arg_type = Int
			default = 1
        "--EVALCONTINUE"
            help = "Conduct evaluation on the final solution concluded as a post process."
            action = :store_true
		"--EVALFILE", "--ef"
			help = "Point to a json file that contains scenarios for evaluation. Search in input folder. If not given, conduct trivia evaluation."
			arg_type = AbstractString
			default = ""
        "--STOCHWRITE"
            help = "Output stochastic file in .json format"
            action = :store_true
        "--DESIGNFILE", "--df"
            help = "Design file (.json) for evaluation. Insert with suffix."
            arg_type = AbstractString
            default = ""
        "--PARAMFILE", "--pf"
            help = "Use input parameters through an json file/parameter generator function"
            arg_type = AbstractString
            default = ""
        "--CGHEURISTIC", "--cg"
            help = "Choose which column generation heuristic to use for sbd_heuristic algorithm"
			arg_type = AbstractString
			default = "improver_heu"
        "--LOCALCONFIG", "--lc"
            help = "Detect configuration file locally inside the repo."
			action = :store_true
        "--HEURISTIC", "--heu"
            help = "Choose which heuristic method for generating adapation plans"
            arg_type = AbstractString
            default = "reactor"
		"--CGMAX"
			help = "Maximum no incumbent out limit on utilizing CG heuristics"
			arg_type = Int
			default = -1
        "--PARAMWRITE"
            help = "Output the current parameter setting in output path in .json format"
            action = :store_true
        "--EVALOBJ", "--eo"
            help = "Evaluation objective function. (| feasibility | slackness |)"
            arg_type = AbstractString
            default = "none"
        "--NAME"
            help = "Output file naming. (Default: 00, you will find design_00.json in output folder.)"
            arg_type = AbstractString
            default = "00"
        "--REPEAT"
            help = "Repeat resampling process with changing random numbers"
            arg_type = Int
            default = 1
        "--CACHE"
            help = "Store/Utilized pre-calculated information. Pre-determined using name code."
            action = :store_true
        "--COSTLambda"
            help = "Tuning parameter [0 to 1]. Cost percentage on expansion. Default is stick with the default value. If set to 1, indicates expensive expansion cost and vise versa"
            arg_type = Float64
            default = -1.0
        "--DEMANDLambda"
            help = "Tuning parameter. Linear demad growth controler. (Default: 0.03)"
			arg_type = Float64
			default = 0.02
        "--SHEDLambda"
            help = "Tuning parameter. Sheeding percentage 5% if 0.95. (Default: 0.95)"
            arg_type = Float64
            default = 0.95
		"--CONGESTLambda"
			help = "Tuning parameter. Adjust thermal limits of all lines. [0-1] CONGESTion ratio (Default: 1.0)"
			arg_type = Float64
			default = 1.0
		"--USESBDNORISK"
			help = "Tuning parameter. Adjust solving method for large no-risk problems"
			arg_type = AbstractString
			default = "config"
		"--ANGLESHIFTLambda"
			help = "Tuning parameter. Angle shift limits. [0-1] of pi. (Default: 0.25)"
			arg_type = Float64
			default = 30.0
		"--DISCOUNTLambda"
			help = "Tuning parameter. Cost discount factor over time. [Default: decreasing discount 1%] (Default: -0.01)"
			arg_type = Float64
			default = -0.01
		"--SEED", "--seed"
			help = "Alternative random seed to replace configuration input"
			arg_type = AbstractString
			default = "config"
        "--SUBSETS"
            help = "Parameter used for resampling from the uncertainty space"
			arg_type = Int
			default = 1
		"--CSV"
			help = "Additional CSV file for assistive information in parameter readings"
			arg_type = AbstractString
			default = ".csv"
		"--DEBUG"
			help = "TURN-ON DEBUGGING mode"
			action = :store_true
		"--OUTPUTPATH", "--out"
			help = "Addditional path director"
			arg_type = AbstractString
			default = ""
    end

    args = parse_args(s)

    # For more flexible usage, override with function inputs
    (haskey(inputDriver, :PROBLEM)) && (args["PROBLEM"] = inputDriver[:PROBLEM])
	(haskey(inputDriver, :MODEL)) && (args["MODEL"] = inputDriver[:MODEL])
	(haskey(inputDriver, :ALG)) && (args["ALG"] = inputDriver[:ALG])
	(haskey(inputDriver, :T)) && (args["T"] = inputDriver[:T])
	(haskey(inputDriver, :S)) && (args["S"] = inputDriver[:S])
	(haskey(inputDriver, :EPS)) && (args["EPS"] = inputDriver[:EPS])
	(haskey(inputDriver, :STOCHMODE)) && (args["STOCHMODE"] = inputDriver[:STOCHMODE])
	(haskey(inputDriver, :STOCHFILE)) && (args["STOCHFILE"] = inputDriver[:STOCHFILE])
	(haskey(inputDriver, :DESIGNFILE)) && (args["DESIGNFILE"] = inputDriver[:DESIGNFILE])
	(haskey(inputDriver, :PARAMFILE)) && (args["PARAMFILE"] = inputDriver[:PARAMFILE])
	(haskey(inputDriver, :EVALOBJ)) && (args["EVALOBJ"] = inputDriver[:EVALOBJ])
	(haskey(inputDriver, :FEATURES)) && (args["FEATURES"] = inputDriver[:FEATURES])
	(haskey(inputDriver, :DEMANDLambda)) && (args["DEMANDLambda"] = inputDriver[:DEMANDLambda])
	(haskey(inputDriver, :SHEDLambda)) && (args["SHEDLambda"] = inputDriver[:SHEDLambda])
	(haskey(inputDriver, :CONGESTLambda)) && (args["CONGESTLambda"] = inputDriver[:CONGESTLambda])
	(haskey(inputDriver, :ANGLESHIFTLambda)) && (args["ANGLESHIFTLambda"] = inputDriver[:ANGLESHIFTLambda])
	(haskey(inputDriver, :DISCOUNTLambda)) && (args["DISCOUNTLambda"] = inputDriver[:DISCOUNTLambda])
	(haskey(inputDriver, :COSTLambda)) && (args["COSTLambda"] = inputDriver[:COSTLambda])

    return args
end

"""
    Additional Parser to parse and analyze user-inputs for driving main function
    Take parsed arguments and default values inherited in the problem
"""
function get_driver_args(args::Dict; kwargs...)

    driver = Dict(kwargs)
    driver[:PROBLEM] = args["PROBLEM"]
    driver[:CACHE] = args["CACHE"]

    (args["T"] > 0) && (driver[:T] = args["T"])
	(args["S"] > 0) && (driver[:S] = args["S"])
    (args["EPS"] != 1.0) && (driver[:eps] = args["EPS"])

    # Stochasticity arguments
    if isempty(args["STOCHFILE"])
        # Check if user indicated any scenario generation method
        (!isempty(args["STOCHMODE"])) && (driver[:STOCHMODE] = args["STOCHMODE"])
        driver[:STOCHFILE] = ""
        driver[:STOCHWRITE] = args["STOCHWRITE"]
    else
        driver[:STOCHMODE] = args["STOCHMODE"]
        driver[:STOCHFILE] = args["STOCHFILE"]
        driver[:STOCHWRITE] = args["STOCHWRITE"]
    end

    isempty(args["PARAMFILE"]) ? driver[:PARAMFILE] = "" : driver[:PARAMFILE] = args["PARAMFILE"]

    driver[:PARAMWRITE] = args["PARAMWRITE"]
    driver[:ALGO] = args["ALG"]
	driver[:REPEAT] = args["REPEAT"]
	driver[:CGHEURISTIC] = args["CGHEURISTIC"]
	driver[:HEURISTIC] = args["HEURISTIC"]
    driver[:EVALDESIGN] = args["DESIGNFILE"]
    driver[:EVALTARGET] = args["EVALOBJ"]
	driver[:EVALCONTINUE] = args["EVALCONTINUE"]
	driver[:EVALFILE] = args["EVALFILE"]
    driver[:NAME] = args["NAME"]
    driver[:PARALLEL] = args["PARALLEL"]
	driver[:WARMSTART] = args["WARMSTART"]
    driver[:SOLVER] = args["SOLVER"]
	driver[:WORKERS] = args["WORKERS"]
    driver[:COSTLambda] = args["COSTLambda"]
    driver[:SHEDLambda] = args["SHEDLambda"]
	driver[:CONGESTLambda] = args["CONGESTLambda"]
	driver[:ANGLESHIFTLambda] = args["ANGLESHIFTLambda"]
	driver[:DISCOUNTLambda] = args["DISCOUNTLambda"]
	driver[:DEMANDLambda] = args["DEMANDLambda"]
	driver[:LEVEL] = args["LEVEL"]
	driver[:SUBSETS] = args["SUBSETS"]
	driver[:CGMAX] = args["CGMAX"]
    driver[:OUTPATH] = args["OUTPUTPATH"]

	if isfile("$(config.INPUTPATH)$(args["CSV"]).csv")
		driver[:CSV] = readtable("$(config.INPUTPATH)$(args["CSV"]).csv")
		info("Found CSV file as support data")
	elseif isfile("$(config.INPUTPATH)$(args["CSV"])")
		driver[:CSV] = readtable("$(config.INPUTPATH)$(args["CSV"])")
		info("Found CSV file as support data")
	else
		driver[:CSV] = nothing
	end

    if args["MODEL"] == "network"
       driver[:MODEL] = network_characteristic
    elseif args["MODEL"] == "capacity"
       driver[:MODEL] = capacity_characteristic
    elseif args["MODEL"] == "dc"
       driver[:MODEL] = dc_characteristic
    elseif args["MODEL"] == "ac"
       error("ERROR|adcc.jl|main()|AC model characteristic not yet implemented.")
    else
       error("ERROR|adcc.jl|main()|Unkown model type. Check input")
    end

    (haskey(args, "FEATURES")) && (driver[:FEATURES] = args["FEATURES"])

	return driver
end
