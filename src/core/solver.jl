# This function returns a empty model with configured solver attached to it
function identify_solver(driver::Dict)

    solverstring = string(driver[:SOLVER])

    if contains(solverstring,"Gurobi")
        return "Gurobi"
    elseif contains(solverstring,"CPLEX")
        return "CPLEX"
    elseif contains(solverstring,"Cbc")
        return "Cbc"
    elseif contains(solverstring,"GLPK")
        return "GLPK"
    else
        error("Unsupported mip solver.")
    end

end

function insert_solver_option(options, kw::Symbol, val::Any)

    for (i,k) in enumerate(options)
        kw in collect(k) && deleteat!(options, i)
    end
    push!(options, (kw, val))

    return
end

function config_solver(model::JuMP.Model, driver::Dict;
                       focus=0,
                       timelimit=-1,
                       showlog=0,
                       mipgap=-1,
                       threads=-1,
                       presolve=1,
                       use_license=false)

    if focus == "feasibility" #1
        focus = 1
    elseif focus == "optimality" #2
        focus = 2
    elseif focus == "bound" #3
        focus = 3
    end

    if timelimit < 0.0
        timelimit = driver[:TIMELIMIT]
    end

    if mipgap < 0.0
        mipgap = driver[:OPTGAP]
    end

    if threads < 0.0
        mipgap = driver[:THREADS]
    end

    solver = driver[:SOLVER]() # Generate a solver
    if identify_solver(driver) == "CPLEX"
        insert_solver_option(solver.options, :CPX_PARAM_TILIM, timelimit)
        insert_solver_option(solver.options, :CPX_PARAM_SCRIND, showlog)
        insert_solver_option(solver.options, :CPX_PARAM_EPGAP, mipgap)
        insert_solver_option(solver.options, :CPX_PARAM_PREIND, presolve)
        insert_solver_option(solver.options, :CPX_PARAM_THREADS, threads)
        insert_solver_option(solver.options, :CPX_PARAM_MIPEMPHASIS, focus)
    elseif identify_solver(driver) == "Gurobi"
        insert_solver_option(solver.options, :TimeLimit, timelimit)
        insert_solver_option(solver.options, :OutputFlag, showlog)
        insert_solver_option(solver.options, :MIPGap, mipgap)
        insert_solver_option(solver.options, :Presolve, presolve)
        insert_solver_option(solver.options, :Threads, threads)
        insert_solver_option(solver.options, :MIPFocus, focus)
        insert_solver_option(solver.options, :SimplexPricing, 1)
        insert_solver_option(solver.options, :Heuristics, 0.001)
    elseif identify_solver(driver) == "Cbc"
        warn("No special support for Cbc Solver options")
    elseif identify_solver(driver) == "GLPK"
        warn("No special support for GLPK Solver options")
    end

    setsolver(model, solver)

    return model
end

function print_iis_gurobi(m::Model, driver::Dict)

	if identify_solver(driver[:SOLVER]) == "Gurobi"

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
