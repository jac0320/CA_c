"""
    Bender's decomposition main driver.
    Iterative method is utilized here.
    There is some issues with the callback method.

    Probelm Structure ::
        Min   cx + eta
        s.t.   Ax = b
        where eta = E(h(x))
            h(x) = min  qy
            s.t.  Wy = h - Tx
"""
function benders(master_variables, master_constraints, master_objective, subprob_variables, subprob_constraints, subprob_objective, characteristic_formulation; kwargs...)

    # ====================== Pre-process ======================= #
    dataPackage = Dict(kwargs)

    master_suite = Dict()
    subprob_suite = Dict()

    master_suite[:variables] = master_variables
    master_suite[:constraints] = master_constraints
    master_suite[:objective] = master_objective

    subprob_suite[:variables] = subprob_variables
    subprob_suite[:constraints] = subprob_constraints
    subprob_suite[:objective] = subprob_objective

    subprob_suite[:characteristic] = characteristic_formulation

    info("Finished setting up the master/subprob suites...")

    # ===================== Initialization ===================== #
    cuts = []
    iter = 0
    S = dataPackage[:stoc].S
    incumbentObj = Inf

    # => Construct master formulation that stores everythign
    # => Initial Master problem
    master = benders_master(master_suite, dataPackage)
    # print(master.model)

    # => Fetch critical matrices from subproblem structure
    # => Think about how this can be better fit into the data structurequ
    GAB_time = @elapsed TMs, hVs, stoMaps = get_analytic_blocks(master, subprob_suite, dataPackage)

    # => Solve master problem (has no eta in objective if eta has no lower bound)
    if getlowerbound(master.vars[:eta]) == -Inf
        master.model = master_suite[:objective](master.model, master.vars, hasEta=false)
    end
    status = solve(master.model)

    # Get master solution in solnType
    masterObj = master.model.objVal
    masterSol = get_primal_solution(master)     # Stores in the solnType
    masterVal = master.model.colVal

    xVector = get_x_vector(master.vars, master.cols) # of symbols

    feaCutCnt = 0
    optCutCnt = 0

    info("Collecting master problem variables (symbolized)...")

    while true
        print("\n")
        info("--- Iteration $iter ---")
        iter += 1

        feaCut = false
        pπh = 0.0
        pπT = zeros(Float64, master.cols)

        for s in 1:S

            # @show dataPackage[:stoc].scenarios[s].data["SS"]
            # Setup Benders Subproblem
            subprob = bender_subprob(subprob_suite,
                                    master.varBuilder,
                                    masterSol,
                                    dataPackage[:stoc].scenarios[s],
                                    dataPackage)

            #@show dataPackage[:stoc].scenarios[s].data

            # Solve subproblem to generate lower bounds => subproblem is lp
            status = solve(subprob.model, suppress_warnings=true)

            # Fetch corresponding transfer matrix
            Tbar = TMs[s]
            hBar = hVs[s]
            hTilde, Ttilde = get_tilde(stoMaps[s], dataPackage[:stoc].scenarios[s], subprob.rows, master.cols)

            if status == :Optimal # Generate Optimality Cut
                print("+")
                if !feaCut
                    # info("Problem feasible. Accumulating optimality cuts...")

                    subprimal = subprob.model.colVal
                    subprimalObj = getobjectivevalue(subprob.model)
                    subprobSol = get_primal_solution(subprob)

                    # print("\n")
                    # info("[Scenario $s] Subprob objective is $subprimalObj")

                    π = subprob.model.linconstrDuals
                    p = dataPackage[:stoc].scenarios[s].chance
                    h = hBar + hTilde
                    T = Tbar + Ttilde

                    validRHS = get_rhs(subprob.model)
                    vectorRHS = h - T * masterVal
                    for v in 1:length(validRHS)
                        @assert validRHS[v] - vectorRHS[v] <= config.TOLERANCE
                    end

                    # ----- Debugging Area ----- #
                    # @show stoMaps[s]
                    # @show hBar
                    # @show hTilde
                    # @show Tbar
                    # @show Ttilde
                    # @show masterSol
                    # @show masterVal
                    # @show subprobSol
                    # @show masterSol
                    # @show (h-T*masterVal)
                    # @show (π'*(h-T*masterVal))[1]
                    # print(subprob.model)
                    # --------------------------- #

                    # Duality Check for Debugging
                    @assert (π'*(h-T*masterVal))[1] - subprimalObj <= 0.1

                    if (s == 1)
                        pπh = p * (π' * h)
                        pπT = p * (π' * T)
                    else
                        pπh += p * (π' * h)
                        pπT += p * (π' * T)
                    end

                end
            elseif status == :Infeasible # Generate Feasibility Cuts
                α = subprob.model.linconstrDuals

                h = hBar + hTilde
                T = Tbar + Ttilde

                validRHS = get_rhs(subprob.model)
                vectorRHS = h - T * masterVal

                for v in 1:length(validRHS)
                    @assert validRHS[v] - vectorRHS[v] <= config.TOLERANCE
                end

                αh = α' * h
                αT = α' * T

                @constraint(master.model, feaCut, dot(vec(xVector),vec(αT)) .>= αh)
                print("^")
                push!(cuts, feaCut)
                feaCut = true
                feaCutCnt += 1
            else
                error("ERROR|bender.jl|add_bender_cut()|Unkown subproblem status.")
            end
        end

        # Master problem should have variable eta available
        if !feaCut   #If no feasibility cuts is added, then add a optimality cut
            @constraint(master.model, optCut, master.model.varDict[:eta] .>= pπh - dot(vec(xVector), vec(pπT)))
            push!(cuts, optCut)
            optCutCnt += 1
            @show "OPT Cut $optCutCnt, $optCut"
        end

        # Optimality Check :: Entering condition need to be re-exmained
        if iter > 1 && optCutCnt >= 1 && !feaCut
            η = master.model.colVal[end]
            ub = pπh - dot(pπT,masterVal)
            # println("Checking optimality LB(eta)[$η] and UB[$ub]...")
            # println(masterSol.primal[:ass])
            # @show "---------------------"
            # println(masterSol.primal[:ah])
            # @show "---------------------"
            # println(masterSol.primal[:f])
            # @show "---------------------"
            # println(masterSol.primal[:pg])
            # @show "---------------------"
            # println(masterSol.primal[:h])
            if η >= ub[1]
                println("Optimality reached! Final objective value = $masterObj")
                return master
            end
        end

        # Modify master objective if no optimality cut has ever been generated
        if getlowerbound(master.vars[:eta]) == -Inf && optCutCnt == 0
            master.model = master_suite[:objective](master.model, master.vars, hasEta=false)
        else
            master.model = master_suite[:objective](master.model, master.vars, hasEta=true)
        end

        # info("Solving master problem ...")
        # setsolver(master.model, GurobiSolver(Presolve=1))
        status = solve(master.model)

        # Get master solution in solnType
        masterObj = master.model.objVal
        masterSol = get_primal_solution(master)     # Stores in the solnType
        masterVal = master.model.colVal
        masterTime = round.(getsolvetime(master.model),2)

        if masterObj < incumbentObj
            incumbentObj = masterObj
            print("+")
        end

        print("\n")
        info("[obj=$masterObj][time=$masterTime] Collecting master solutions ...")
    end
end

"""
	This master formulation is different from the regular formulation constructor format. Pay sepecial attention.
"""
function benders_master(master_suite::Dict, dataPackage::Dict; kwargs...)

    options = Dict(kwargs)

    m = init_model_solver()

    # Construct variables
    m, vars = master_suite[:variables](m, package=dataPackage)
    m = master_suite[:constraints](m, vars, package=dataPackage)
    m = master_suite[:objective](m, vars, hasEta=true)

    master = oneProblem()
    master.name = "BendersMaster"
    master.T = dataPackage[:stoc].T
    master.S = dataPackage[:stoc].S
    master.stage = 1
    master.cols = m.numCols
    master.rows = MathProgBase.numlinconstr(m)
    master.vars = vars
    master.param = dataPackage[:param]
    master.model = m
    master.samples = dataPackage[:stoc]
    master.builder = nothing
    master.varBuilder = master_suite[:variables]
    master.consBuilder = master_suite[:constraints]
    master.objBuilder = master_suite[:objective]
    master.status = :Unsolved
    master.objective = Inf

	return master
end

"""
	Description about bender's subprob problem goes here.
"""
function bender_subprob(subprob_suite, master_variables, masterX::solnType, scenario::scenarioType, dataPackage::Dict, analytic::Bool=false; kwargs...)

    options = Dict(kwargs)

    m = init_model_solver()
    scenIdx = scenario.ind

    # Here, the sequence matters
    if analytic
        m, Xvars = master_variables(m, package=dataPackage)
        masterColCnt = m.numCols    #Eta included
    else
        Xvars = masterX.primal
    end

    m, Yvars  = subprob_suite[:variables](m, package=dataPackage)
    m, stoMap = subprob_suite[:constraints](m, Yvars, Xvars, scenario, package=dataPackage)
    m, stoMap = subprob_suite[:objective](m, Yvars, scenario, stoMap, package=dataPackage)

    if analytic
        TMatrix = benders_get_Tbar(m, Xvars, masterColCnt)     #Also TBar
        hBar = benders_get_hBar(m)
        return TMatrix, hBar, stoMap
    end

	subprob = oneProblem()
    subprob.name = string("BendersSubprob_",scenIdx)
    subprob.T = dataPackage[:stoc].T
    subprob.S = 1
    subprob.eps = 0.0
    subprob.stage = 2
    subprob.cols = m.numCols
    subprob.rows = MathProgBase.numlinconstr(m) # Potentially can have SOC
    subprob.vars = Yvars
    subprob.param = dataPackage[:param]
    subprob.model = m
    subprob.samples = scenario
    subprob.builder = nothing
    subprob.varBuilder = subprob_suite[:variables]
    subprob.consBuilder = subprob_suite[:constraints]
    subprob.objBuilder = subprob_suite[:objective]
    subprob.status = :Unsolved
    subprob.objective = Inf

	return subprob

end

"""
    Input subprob must be a model that has master problem symbols.
"""
function benders_get_Tbar(spm::JuMP.Model, Xvars::Dict, masterColCnt::Int; kwargs...)

    I = Int[]
    J = Int[]
    V = Float64[]

    subRCnt = MathProgBase.numlinconstr(spm)

    # TODO: this needs to go to parallel
    for (nrow, con) in enumerate(spm.linconstr)
        aff = con.terms
        for (var, id) in zip(reverse(aff.vars), length(aff.vars):-1:1)
            isMasterVar = false
            for varCluster in keys(Xvars)
                if var in Xvars[varCluster]
                    isMasterVar = true
                end
            end
            if isMasterVar
                # Tbar = add2SparseMatrix(Tbar, nrow, var.col, aff.coeffs[id])
                push!(I,nrow)
                push!(J,var.col)
                push!(V,aff.coeffs[id])
            end
        end
    end

    return sparse(I,J,V,subRCnt,masterColCnt)
end

function benders_get_hBar(spm::JuMP.Model)

    # hBar = sparseVector()
    I = Int[]
    V = Float64[]

    subRCnt = MathProgBase.numlinconstr(spm)

    for (nrow, con) in enumerate(spm.linconstr)
        # Be careful, there shouldn't be any two sided constraints
        # Be careful, these only allows one random value at one constraint at a time ???
        if (con.lb == 0 && con.ub == Inf) || (con.lb == -Inf && con.ub == 0) || (con.lb == 0 && con.ub == 0)
            skip = 0    # Zero bound case
        elseif (con.lb == -Inf) # >= constraint with non-zero rhs
            # hBar = add2SparseVector(hBar, nrow, con.ub)
            push!(I, nrow)
            push!(V, con.ub)
        elseif (con.ub == Inf)  # >= constraint with non-zero rhs
            # hBar = add2SparseVector(hBar, nrow, con.lb)
            push!(I, nrow)
            push!(V, con.lb)
        else                    # Equality constraint with non-zero rhs
            # hBar = add2SparseVector(hBar, nrow, con.lb)
            push!(I, nrow)
            push!(V, con.lb)
        end
    end

    return sparsevec(I, V, subRCnt)
end

"""
    This function takes in an realization of stochasticity and problem's stochasticity map and generates the tilde part of vector h and transfer matrix T for cut generation.
"""
function get_tilde(stoMap::Dict, scenario::scenarioType, subRCnt::Int, masterCCnt::Int)

    # The tilde part for h
    # hTilde = sparseVector()
    Ih = Int[]
    Vh = Float64[]

    stoMap_h = stoMap[:h]
    for stoItem in stoMap_h
        stoName  = stoItem[1][1]
        stoIdx   = Base.CartesianIndex(tuple(stoItem[1][2:1:length(stoItem[1])]...))
        stoValue = scenario.data[stoName][stoIdx]
        rowIdx   = stoItem[2]
        coeffs   = stoItem[3]
        # hTilde = add2SparseVector(hTilde, rowIdx, coeffs*stoValue)
        push!(Ih,rowIdx)
        push!(Vh,coeffs*stoValue)
    end

    # Ttilde = sparseMatrix()
    IT = Int[]
    JT = Int[]
    VT = Float64[]
    stoMap_T = stoMap[:T]
    for stoItem in stoMap_T
        stoName  = stoItem[1][1]
        stoIdx   = Base.CartesianIndex(tuple(stoItem[1][2:1:length(stoItem[1])]...))
        stoValue = scenario.data[stoName][stoIdx]
        rowIdx   = stoItem[2].idx
        colIdx   = stoItem[3].col
        coeffs   = stoItem[4]
        # Ttilde = add2SparseMatrix(Ttilde, rowIdx, colIdx, coeffs*stoValue)
        push!(IT, rowIdx)
        push!(JT, colIdx)
        push!(VT, coeffs*stoValue)
    end

    # return hTilde, Ttilde
    return sparsevec(Ih, Vh, subRCnt), sparse(IT, JT, VT, subRCnt, masterCCnt)
end

function get_x_vector(Xvars::Dict, Xcnt::Int)

    xVector = Array{JuMP.Variable}(Xcnt)

    for varClusterKey in keys(Xvars)
        for item in Xvars[varClusterKey]
            xVector[item.col] = item
        end
    end

    return xVector

end

"""
    This function retrive matrices h, T, q
"""
function get_analytic_blocks(master::oneProblem, subprob_suite::Dict, dataPackage::Dict; kwargs...)

    info("Obtaining analytical matrices/vectors...")
    S = dataPackage[:stoc].S

    @show S

    TMs     = Array{SparseMatrixCSC}(S)
    hVs     = Array{SparseVector}(S)
    stoMaps = Array{Dict}(S)

    #
    # reducer = Array{Tuple}(S)
    # cStoc   = copy_null_stocType(dataPackage[:stoc])
    #
    #
    # reducer = pmap((a1,a2,a3,a4,a5,a6)->bender_subprob(a1,a2,a3,a4,a5,a6),
    #                         [subprob_suite for s in 1:S],
    #                         [master.varBuilder for s in 1:S],
    #                         [solnType() for s in 1:S],
    #                         [cStoc.scenarios[s] for s in 1:S],
    #                         [dataPackage for s in 1:S],
    #                         [true for s in 1:S])
    #
    # # Remap to accurate locations
    # for s in 1:S
    #     TMs[s] = reducer[s][1]
    #     hVs[s] = reducer[s][2]
    #     stoMaps[s] = reducer[s][3]
    # end

    # TMs     = Dict()
    # hVs     = Dict()
    # stoMaps = Dict()

    cStoc = copy_null_stocType(dataPackage[:stoc])

    for s in 1:dataPackage[:stoc].S
        TM, hV, stoMap = bender_subprob(subprob_suite,
                                        master.varBuilder,
                                        solnType(),
                                        cStoc.scenarios[s],
                                        dataPackage,
                                        true)
        TMs[s] = TM
        hVs[s] = hV
        stoMaps[s] = stoMap
    end

    return TMs, hVs, stoMaps
end
