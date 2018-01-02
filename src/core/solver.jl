# This function returns a empty model with configured solver attached to it
function init_model_solver()

    if config.SOLVER == "CPLEX" || config.SOLVER == "Cplex"
        m = Model(solver = CplexSolver())
        solver_config(m)
    elseif config.SOLVER == "GUROBI" || config.SOLVER == "Gurobi"
        m = Model(solver = GurobiSolver())
        solver_config(m)
    elseif config.SOLVER == "Cbc"
        m = Model(solver = CbcSolver())
    else
		error("ERROR::solver.jl::init_model_solver()::Unsupport Solver... Change configuration file or give command line argument --SOLVER");
    end

    return m
end

function solver_config(model::JuMP.Model; kwargs...)

    options = Dict(kwargs)

    if haskey(options, :focus)
        if options[:focus] == "feasibility" #1
            focus = 1
        elseif options[:focus] == "optimality" #2
            focus = 2
        elseif options[:focus] == "bound" #3
            focus = 3
        else
            info("Unkown focus option. Setting to 0 - balanced.")
            focus = 0
        end
    else
        focus = 0
    end

    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = config.TIMELIMIT
    haskey(options, :showlog) ? showlog = options[:showlog] : showlog = 0
    haskey(options, :mipgap) ? mipgap = options[:mipgap] : mipgap = config.OPTGAP
    haskey(options, :threads) ? threads = options[:threads] : threads = config.THREADS
    haskey(options, :presolve) ? presolve = options[:presolve] : presolve = 1
    haskey(options, :license) ? use_license = true : use_license = false

    if config.SOLVER == "CPLEX" || config.SOLVER == "Cplex"
        setsolver(model, CplexSolver(CPX_PARAM_TILIM=timelimit,
                                       CPX_PARAM_SCRIND=showlog,
                                       CPX_PARAM_EPGAP=mipgap,
                                       CPX_PARAM_PREIND=presolve,
                                       CPX_PARAM_THREADS=threads,
                                       CPX_PARAM_MIPEMPHASIS=focus));
    elseif config.SOLVER == "GUROBI" || config.SOLVER == "Gurobi"
        if use_license && options[:license] != 1
            setsolver(model, GurobiSolver(options[:license],
                                        OutputFlag=showlog,
                                        TimeLimit=timelimit,
                                        MIPGap=mipgap,
                                        Presolve=presolve,
                                        Threads=threads,
                                        MIPFocus=focus,
                                        SimplexPricing=1,Heuristics=0.001));
        else
            setsolver(model, GurobiSolver(OutputFlag=showlog,
                                        TimeLimit=timelimit,
                                        MIPGap=mipgap,
                                        Presolve=presolve,
                                        Threads=threads,
                                        MIPFocus=focus,
                                        SimplexPricing=1,Heuristics=0.001));
        end
    elseif config.SOLVER == "Cbc"
        info("Cbc solver doesn't support MIP focus and Presolve")
        setsolver(model, CbcSolver(seconds=timelimit,
                                    logLevel=showlog,
                                    ratioGap=mipgap,
                                    threads=threads));
    else
        error("ERROR::solver.jl::init_model_solver()::Unsupport Solver... Change configuration file");
    end

    return model
end


#This subroutine writes down the problem in LP/MPS form
function writeProblem(m::JuMP.Model, filename::AbstractString="", fileType::AbstractString="")
	if isempty(filename) == true
		print(m) #Print model (deafult LP) to screen
	else
		if fileType == "LP" || fileType == ""
			writeLP(m, filename);
		elseif fileType == "MPS"
			writeMPS(m, filename);
		else
			print("Unkown problem type. Please check input");
		end
	end
end

function solver_lower_bound(m::JuMP.Model)

    lb = 0.0
    offset = 0.0

	offset = getobjective(m).aff.constant
    lb = getobjbound(m) + offset

	return lb
end

function print_iis_gurobi(m::Model)

	if config.SOLVER == "Gurobi"

	    grb = MathProgBase.getrawsolver(internalmodel(m))
	    Gurobi.computeIIS(grb)
	    numconstr = Gurobi.num_constrs(grb)
	    numvar = Gurobi.num_vars(grb)

	    iisconstr = Gurobi.get_intattrarray(grb, "IISConstr", 1, numconstr)
	    iislb = Gurobi.get_intattrarray(grb, "IISLB", 1, numvar)
	    iisub = Gurobi.get_intattrarray(grb, "IISUB", 1, numvar)

	    info("Irreducible Inconsistent Subsystem (IIS)")
	    info("Variable bounds:")
	    for i in 1:numvar
	        v = Variable(m, i)
	        if iislb[i] != 0 && iisub[i] != 0
	            println(getlowerbound(v), " <= ", getname(v), " <= ", getupperbound(v))
	        elseif iislb[i] != 0
	            println(getname(v), " >= ", getlowerbound(v))
	        elseif iisub[i] != 0
	            println(getname(v), " <= ", getupperbound(v))
	        end
	    end

	    info("Constraints:")
	    for i in 1:numconstr
	        if iisconstr[i] != 0
	            println(m.linconstr[i])
	        end
	    end

	else
		info("Current solver doesn't support IIS.")
	end

end

function setup_envs()

    if config.SOLVER == "Gurobi"
    	config.ENVS = Gurobi.Env()
    else
        config.ENVS = 1
    end

end
