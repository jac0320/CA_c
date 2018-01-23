function get_scenarios(extras::Dict)

    srand(config.RUNSEED); #Set random seed
    if !isempty(extras[:STOCHFILE])
        scenarios = read_stoch_file(extras[:STOCHFILE], extras)
        extras[:S] = scenarios.S
        return scenarios
    else
        error("No scenario file indicated.")
    end
end

function read_stoch_file(filepath::String, exargs::Dict)

    info("Stochastic file input :: $filepath")
    if isfile(filepath)
        sinput = JSON.parsefile(filepath)
    elseif isfile(joinpath(config.INPUTPATH, exargs[:PROBLEM],filepath))
        sinput = JSON.parsefile(joinpath(config.INPUTPATH, exargs[:PROBLEM], filepath))
    elseif isfile(joinpath(config.INPUTPATH, filepath))
        sinput = JSON.parsefile(joinpath(config.INPUTPATH,filepath))
    else
        error("ERROR|stoch.jl|read_stoch_file()|Undetected scenario file. Check your input arguments.")
    end

    # Validatea file intactness :: always require scenarios to provide both SLR & SS
    !haskey(sinput, "SL") && error("No SLR indicated in input file.")
    !haskey(sinput, "SS") && error("No SLR indicated in input file.")

    @assert length(sinput["SL"]) >= 1
    @assert length(sinput["SL"]) == length(sinput["SS"])
    @assert length(sinput["SL"]["1"]) == length(sinput["SS"]["1"])

    S = length(sinput["SL"])
    T = length(sinput["SL"]["1"])
    B = length(sinput["SS"]["1"][1])

    stoc = stocType(S, T, B)

    userT = exargs[:T]
    if userT > stoc.T
        error("Trying to run model with high resolution failed. Scenario files is with less.")
    elseif userT == stoc.T
        info("User model resolution matchs scenario inputs. No additional action needed.")
        userSteps = [1:stoc.T;]
    else
        info("User model resolution coraser. With fixed horizon... resampling by corasen resolutions.")
        userSteps = []
        if stoc.T % userT == 0
            userStepSize = Int(stoc.T/userT)
            for t in 1:userT
                push!(userSteps, t * userStepSize)
            end
            info("Scenario file input resolution $(stoc.T). Request resolution $(userT). ")
            info("Extracting data by steps $(userSteps)")
            info("Overriding stoc.T ")
            stoc.T = userT
        else
            error("Time resolution corasion failed with mismatched data.")
        end
    end

    stoc.scenarios = Array{scenarioType}(stoc.S)

    for i in 1:stoc.S
        s = scenarioType(i,Dict("SL"=>zeros(Float64,T),"SS"=>zeros(Float64,B,T)),1/stoc.S)
        j = 1
        for t in userSteps
            s.data["SL"][j] = sinput["SL"][string(i)][t]
            for b in 1:stoc.B
                s.data["SS"][b,j] = sinput["SS"][string(i)][t][b] # Dimension resize
            end
            j += 1
        end
        stoc.scenarios[i] = s
    end

    return stoc
end

function create_null_samples(B::Int, T::Int)

    samples = stocType()
    samples.S = 1
    samples.T = T
    samples.B = 1

    sl = zeros(Float64,T)
    ss = zeros(Float64,B,T)

    samples.scenarios = []
    push!(samples.scenarios, scenarioType())

    samples.scenarios[1].ind = 1
    samples.scenarios[1].data = Dict("SL"=>sl, "SS"=>ss)
    samples.scenarios[1].chance = 1.0
    samples.scenarios[1].pool = []

    return samples
end

function summary_scenarios(stoc::stocType, param::Dict)

    info("[STOCH] Total Scenario $(stoc.S)")
    info("[STOCH] Total Time Step Count $(stoc.T)")
    info("[STOCH] Data shape => SLR $(size(stoc.scenarios[1].data["SL"]))")
    info("[STOCH] Data shape => SS $(size(stoc.scenarios[1].data["SS"]))")
    info("[STOCH] Minimum Bus Elevation => $(minimum(param[:Ele]))")
    max_slr = 0.0
    for s in 1:stoc.S
        max_slr = max(max_slr, maximum(stoc.scenarios[s].data["SL"]))
    end
    info("[STOCH] MAX-SLR => $(max_slr)")
    info("[STOCH] BUS under MAX-SLR => $(length(param[:Ele][param[:Ele] .<= max_slr]))")
    for s in 1:stoc.S
        info("[STOCH][S=$s] SLR $(stoc.scenarios[s].data["SL"])")
        info("[STOCH][S=$s] BUS under SS ")
        for t in 1:stoc.T
            ss_bus_cnt = 0
            for b in 1:param[:B]
                if param[:Ele][b] < stoc.scenarios[s].data["SS"][b,t]
                    ss_bus_cnt += 1
                    info("BUS $(b) surged by $(round(stoc.scenarios[s].data["SS"][b,t]-param[:Ele][b],2))m without protection. (LOAD = $(param[:Pd][b]) | CAP = $(param[:PgUB][b]))")
                end
            end
            println("T$(t)=>$(ss_bus_cnt); ")
        end
        print("\n")
    end

    return
end

function write_stocType_json(stoc::stocType, filename::AbstractString)
    stocDict = Dict()
    stocDict["SL"] = Dict()
    stocDict["SS"] = Dict()

    for s in 1:stoc.S
        stocDict["SL"][s] = stoc.scenarios[s].data["SL"]
        stocDict["SS"][s] = stoc.scenarios[s].data["SS"]
    end

    wf = open(string(config.OUTPUTPATH,"/",filename), "w")
    write(wf, JSON.json(stocDict))
    close(wf)

    return
end

function subsetting_stocType(stoc::stocType, subset, exargs::Dict)

    # Two main mode :: 1->subset is int; 2->set/array
    if isa(subset, Int)
        if subset < stoc.S
            selected = sample(1:stoc.S,subset,replace=false)
        elseif subset == stoc.S
            return stoc #Nothing need to be done
        else
            error("ERROR|stoch.jl|subsetting_stocType()|Cardinality of the subset cannot be more than original.")
        end
    else
        selected = subset
    end

    # Construct new stocType structure
    substoc = stocType()
    substoc.S = length(selected)
    substoc.T = length(exargs[:T])
    substoc.scenarios = []

    for s in selected
        oneScenario = scenarioType()
        oneScenario.data = stoc.scenarios[s].data
        oneScenario.chance = 1/length(selected)
        oneScenario.pool = []
        push!(substoc.scenarios, oneScenario)
    end

    return substoc
end

function append_scenarios(stoc::stocType, exstoc::stocType)

    # Need to conduct some assertion
    @assert stoc.B == exstoc.B
    @assert stoc.T == exstoc.T

    # Perform the joint
    for s in exstoc.scenarios
        push!(stoc.scenarios, s)
        stoc.S += 1
    end

    return stoc
end

function get_percentile_SS(stoc::stocType, perc::Int, select=[])

    if isempty(select)
        select = collect(1:stoc.S)
    end

    allSS = []
    for s in select
        append!(allSS, stoc.scenarios[s].data["SS"])
    end

    idx = Int(round.(perc/100*length(allSS)))
    res = sort(allSS)[idx]

    return res
end

function compose_scenarios(SSMode::AbstractString, SLRFilePath::AbstractString, driver::Dict)

    B = driver[:B]
    T = driver[:T]
    driver[:STOCHMODE] = "file"

    # Read partial scenarios from a file, SS will be left blank here
    SLRStoc, driver = create_samples(B, driver[:S], T, driver, filepath = SLRFilePath)

    # Now generate paritial scenarios internally
    driver[:STOCHMODE] = SSMode
    SSStoc, driver = create_samples(B, driver[:S], T, driver)

    # Now join two scenarios sets
    for s in 1:driver[:S]
        for b in 1:B
            for t in 1:T
                SLRStoc.scenarios[s].data["SS"][b,t] = SLRStoc.scenarios[s].data["SS"][t] + SSStoc.scenarios[s].data["SS"][b,t]
            end
        end
    end

    return SLRStoc
end

function copy_null_stocType(stoc::stocType)

    S = stoc.S

    newStoc = stocType()
    newStoc.S = stoc.S
    newStoc.T = stoc.T
    newStoc.B = stoc.B
    newStoc.scenarios = []
    for s in 1:S
        oneScenario = scenarioType()
        oneScenario.ind = s
        oneScenario.chance = 1/S
        oneScenario.data = copy(stoc.scenarios[s].data)
        for dataKeys in keys(oneScenario.data)
            oneScenario.data[dataKeys] -= oneScenario.data[dataKeys]
        end
        oneScenario.pool = []
        push!(newStoc.scenarios, oneScenario)
    end

    return newStoc
end
