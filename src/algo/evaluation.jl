function evaluation(power::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

    options = Dict(kwargs)

    outfname = string(config.OUTPUTPATH, exargs[:OUTPATH], "eval_", exargs[:NAME], ".out")
    outf = open(outfname, "w")
    evaltodo = []
    evalnames = []
    if isfile(exargs[:EVALDESIGN])
        write(outf, "[EVALUATING] Detecting a design file...")
        push!(evalnames, exargs[:EVALDESIGN])
        push!(evaltodo, parse_design(exargs[:EVALDESIGN], param))
    elseif isfile(string(config.OUTPUTPATH, exargs[:EVALDESIGN]))
        write(outf, "[EVALUATING] Detecting a design file...")
        push!(evalnames, exargs[:EVALDESIGN])
        push!(evaltodo, parse_design(string(config.OUTPUTPATH, exargs[:EVALDESIGN]), param))
    elseif isdir(exargs[:EVALDESIGN])
        write(outf, "[EVALUATING] Detecting a dir with design json")
        filelist = glob("*.json", string(string(exargs[:EVALDESIGN])))
        for f in filelist
            write(outf, "[EVALUATING] Reading design file $(f)\n")
            push!(evalnames, f)
            push!(evaltodo, parse_design(f, param))
        end
    else
        error("Could not locate design .json file")
    end

    allSubprobs = Array{oneProblem}(stoc.S)
    write(outf, "[EVALUATING] Warming start the testing models...")
    warm_start = time()
    for s = 1:stoc.S
        allSubprobs[s] = oneProblem()
        allSubprobs[s] = sbd_base_formulation(power, param, stoc)
        allSubprobs[s] = attach_scenario(allSubprobs[s], stoc, [s], exargs[:MODEL], 0.0, exargs, subprobType="free")
        @objective(allSubprobs[s].model, Min, 0.0)
    end
    warm_cpu_time = time() - warm_start
    write(outf, "[EVALUATING] Finished warm starting (took: $(warm_cpu_time)s)\n")

    for d in 1:length(evaltodo)
        perform_one_eval(power, param, stoc, evaltodo[d], evalnames[d], exargs, outf, allSubprobs)
    end

    close(outf)
    return
end


function perform_one_eval(power, param, stoc, design, designName, exargs, outf, builtModels=[])

    if exargs[:EVALTARGET] == "none"
        targetObj = "feasibility"
    else
        targetObj = exargs[:EVALTARGET]
    end

    write(outf, "[EVALUATING] Target design file $designName\n")
    write(outf, "[EVALUATION] UB = $(get_design_cost(design, param))\n")
    write(outf, "[EVALUATION] TARGET Objective = $(targetObj)\n")
    write(outf, "[EVALUATION] Problem Mode -> $(exargs[:MODEL])\n")

    objPool = []
    for i in 1:stoc.S
        evalProblem = eval_formulation(power, param, stoc, [i], design, exargs, targetObj, builtModels[i])
        solver_config(evalProblem.model,MIPFocus="feasibility", showlog=0, presolve=1)
        status = solve(evalProblem.model, suppress_warnings=true)
        if targetObj == "feasibility"
            if status == :Optimal
                evalObj = 1
            else
                evalObj = 0
                println("Infeasible : $status $i")
                # print_iis_gurobi(evalProblem.model)
            end
        else
            if status == :Infeasible
                print_iis_gurobi(evalProblem.model)
            else
                evalObj = getobjectivevalue(evalProblem.model)
            end
		end
        push!(objPool, evalObj)
    end

    minObj, minObjIdx = findmin(objPool)
    maxObj, maxObjIdx = findmax(objPool)
    aveObj = mean(objPool)
    sumObj = sum(objPool)

    # Small summary of the isolation stage
    write(outf, "[EVALUATION] Results | total = $sumObj | min = $minObj | average = $aveObj | max = $maxObj\n")

    return
end
