function get_scenarios(exargs::Dict)

    srand(config.RUNSEED); #Set random seed
    if !isempty(exargs[:STOCHFILE])
        scenarios = read_stoch_file(exargs[:STOCHFILE], exargs)
        exargs[:S] = scenarios.S
    else
        error("No scenario file indicated.")
    end

    return reprocess_scenarios(scenarios, exargs)
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

function reprocess_scenarios(stoc::stocType, exargs::Dict)

    exargs[:STOCHMODE] == nothing && return stoc

    if exargs[:STOCHMODE] == "max"
        return downscale_scenarios(stoc, 100)
    elseif exargs[:STOCHMODE] == "min"
        return downscale_scenarios(stoc, 1)
    elseif exargs[:STOCHMODE] == "median"
        return downscale_scenarios(stoc, 50)
    elseif exargs[:STOCHMODE] == "ave"
        return downscale_scenarios(stoc, 0.0)
    elseif ismatch(r"\d+-dev", exargs[:STOCHMODE])
        deviation = Float(split(exargs[:STOCHMODE], "-")[1])
        return downscale_scenarios(stoc, deviation)
    elseif ismatch(r"\d+-perc", exargs[:STOCHMODE])
        perc = Int(split(exargs[:STOCHMODE], "-")[1])
        return downscale_scenarios(stoc, perc)
    else
        error("Unknown stochastic mode.")
    end
end

function downscale_scenarios(stoc::stocType, deviation::Float64)

    s = stocType(1, stoc.T, stoc.B)
    s.scenarios[1] = scenarioType(1, Dict("SL"=>zeros(Float64,stoc.T),"SS"=>zeros(Float64,stoc.B,stoc.T)), 1.0)

    slr_sort = zeros(Float64, stoc.S)
    ss_sort = zeros(Float64, stoc.S)

    for t in 1:stoc.T
        for s in 1:stoc.S
            slr_sort[s] = stoc.scenarios[s].data["SL"][t]
        end
        slr_std = std(slr_sort)
        slr_mean = mean(slr_sort)
        s.scenarios.data["SL"][t] = slr_mean + deviation*slr_std

        for b in 1:B
            for s in 1:stoc.S
                ss_sort = stoc.scenarios[s].data["SS"][b,t]
            end
            s.scenario.data["SS"][b,t] = select(sstemp_bnew, Int(rank))
            ss_std = std(ss_sort)
            ss_mean = mean(ss_sort)
            s.scenario.data["SS"][b,t] = ss_mean + deviation*ss_std
        end
    end

    return s
end

function downscale_scenarios(stoc::stocType, perc::Int)

    rank = max(floor((perc/100) * stoc.S), 1)

    s = stocType(1, stoc.T, stoc.B)
    s.scenarios[1] = scenarioType(1, Dict("SL"=>zeros(Float64,stoc.T),"SS"=>zeros(Float64,stoc.B,stoc.T)), 1.0)

    slr_sort = zeros(Float64, stoc.S)
    ss_sort = zeros(Float64, stoc.S)

    for t in 1:stoc.T
        for s in 1:stoc.S
            slr_sort[s] = stoc.scenarios[s].data["SL"][t]
        end
        s.scenarios.data["SL"][t] = select(sltempnew, rank)

        for b in 1:B
            for s in 1:stoc.S
                ss_sort = stoc.scenarios[s].data["SS"][b,t]
            end
            s.scenario.data["SS"][b,t] = select(sstemp_bnew, Int(rank))
        end
    end

    return s
end

function null_scenario_stoc(B::Int, T::Int)

    s = stocType(1, T, B)
    s.scenarios[1] = scenarioType(1, Dict("SL"=>zeros(Float64,stoc.T),"SS"=>zeros(Float64,stoc.B,stoc.T)), 1.0)

    return s
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

    outputDict = Dict()

    outputDict["SL"] = Dict()
    outputDict["SS"] = Dict()

    for s in 1:stoc.S
        outputDict["SL"][s] = stoc.scenarios[s].data["SL"]
        outputDict["SS"][s] = stoc.scenarios[s].data["SS"]
    end

    f = open(joinpath(config.OUTPUTPATH, filename), "w")
    JSON.print(f, outputDict)
    close(f)

    return
end

function subset_scenarios(stoc::stocType, subset::Vector, exargs::Dict)

    N = length(subset)
    s = stocType(N, stoc.T, stoc.B)
    i = 1
    for s in subset
        scenarioType(i, Dict("SL"=>copy(stoc.scenarios[s].data["SL"]),
                             "SS"=>copy(stoc.scenarios[s].data["SS"])), 1/N)
        i += 1
    end

    return s
end

subset_scenarios(stoc::stocType, subset::Set, exargs::Dict) = subset_scenarios(stoc, subset, exargs)

function subset_scenarios(stoc::stocType, subset::Int, exargs::Dict)

    subset > stoc.S && error("subsetting more than original set")
    selected = randperm(stoc.S)[1:subset]

    return subset_scenarios(stoc, selected, exargs)
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

function copy_null_stocType(s::stocType)

    ns = stocType(s.S, s.T, s.B)
    for i in 1:s.S
        ns.scenarios[i] = scenarioType(i, Dict("SL"=>copy(s.scenarios[i].data["SL"]),
                                               "SS"=>copy(s.scenarios[i].data["SS"])), 1/s.S)
    end

    return ns
end
