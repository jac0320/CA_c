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

    println("Stochastic file input :: $filepath")
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
        println("User model resolution matchs scenario inputs. No additional action needed.")
        userSteps = [1:stoc.T;]
    else
        println("User model resolution coraser. With fixed horizon... resampling by corasen resolutions.")
        userSteps = []
        if stoc.T % userT == 0
            userStepSize = Int(stoc.T/userT)
            for t in 1:userT
                push!(userSteps, t * userStepSize)
            end
            println("Scenario file input resolution $(stoc.T). Request resolution $(userT). ")
            println("Extracting data by steps $(userSteps)")
            println("Overriding stoc.T ")
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

# """
#     Input -> B::Int | S::Int | T::Int | sampleMode::Array | ...
#     This is main simulator shell that manage randomness in problems.
#     Output -> stoc::stocType
# """
# function create_samples(B::Int, S::Int, T::Int, exargs; kwargs...)
#
#     extras = Dict(kwargs)
#
#     sampleMode = exargs[:STOCHMODE]
#
#     srand(config.RUNSEED); #Set random seed
#
#     samples = stocType();
#     samples.S = S  # This will be correct if reading from a file
#     samples.T = T
#     samples.B = B
#     samples.sbdColumns = []
#
#     # Driver Arguments
# 	mode = exargs[:STOCHMODE]
#
#     function one_null_scenario()
#         oneScenario = scenarioType()
#         oneScenario.data = Dict("SL"=>zeros(Float64, T), "SS"=>zeros(Float64, B, T))
#         oneScenario.chance = 0.0
#         oneScenario.pool = []
#         return oneScenario
#     end
#
#     function generate_zero_sample(stoc::stocType)
#         stoc.scenarios = []
#         for s in 1:S
#             push!(stoc.scenarios, scenarioType())
#             sl = zeros(Float64,T)
#             ss = zeros(Float64,B,T)
#             push!(stoc.scenarios, scenarioType())
#             stoc.scenarios[s].ind = s
#             stoc.scenarios[s].data = Dict("SL"=>sl, "SS"=>ss)
#             stoc.scenarios[s].chance = 1/S
#             stoc.scenarios[s].pool = []
#         end
#         return stoc
#     end
#
#     function generate_random_sample(stoc::stocType, SLINIT::Float64, serverity::Float64)
#         stoc.scenarios = []
#         for s in 1:S
#             sl = zeros(Float64,T);
#             ss = zeros(Float64,B,T);
#             push!(stoc.scenarios, scenarioType())
#             # Scenario data generation
#             sl[1] = SLINIT + rand()*2
#             for t in 2:T
# 				sl[t] = sl[t-1] + rand()*2;
# 			end
#             for t in 1:T
# 				ss[:,t] = sl[t] + rand(B)*serverity;
# 			end
#             stoc.scenarios[s].ind = s
#             stoc.scenarios[s].data = Dict("SL"=>sl, "SS"=>ss)
#             stoc.scenarios[s].chance = 1/S
#             stoc.scenarios[s].pool = []
#         end
#         return stoc
#     end
#
#     function generate_evolving_sample(stoc::stocType, SLINIT::Float64, serverity_base::Float64)
#         stoc.scenarios = []
#         for s in 1:S
#             sl = zeros(Float64,T)
#             ss = zeros(Float64,B,T)
#             push!(stoc.scenarios, scenarioType())
#             # Scenario data generation
#             sl[1] = SLINIT + rand();
#             for t in 2:T
#                 sl[t] = sl[t-1] + rand() * t * 0.2;
#             end
#             for t in 1:T
#                 randTemp = rand(Exponential(t*serverity_base/2))
#                 for b in 1:B
#                     ss[b,t] = sl[t] + max(0, randTemp);
#                 end
#             end
#             stoc.scenarios[s].ind = s
#             stoc.scenarios[s].data = Dict("SL"=>sl, "SS"=>ss)
#             stoc.scenarios[s].chance = 1/S
#             stoc.scenarios[s].pool = []
#         end
#         return stoc
#     end
#
#     function generate_tuning_sample(stoc::stocType)
#         stoc.scenarios = []
#         for s in 1:S
#             sl = zeros(Float64,T)
#             ss = zeros(Float64,B,T)
#             push!(stoc.scenarios, scenarioType())
#             # Scenario data generation
#             sl[1] = 0.0;
#             for t in 2:T
#                 sl[t] = sl[t-1] + rand() * t * 0.2;
#             end
#             for t in 1:T
#                 randTemp = abs(rand(Normal(2.0, 0.2*t)))
#                 for b in 1:B
#                     ss[b,t] = max(0, randTemp);
#                 end
#             end
#             stoc.scenarios[s].ind = s
#             stoc.scenarios[s].data = Dict("SL"=>sl, "SS"=>ss)
#             stoc.scenarios[s].chance = 1/S
#             stoc.scenarios[s].pool = []
#         end
#         return stoc
#     end
#
#     function calculate_mean_sample(stoc::stocType)
#         meanScenario = scenarioType()
#         meanScenario.ind = 1
#         meanScenario.data = Dict("SL"=>zeros(Float64, T), "SS"=>zeros(Float64, B, T))
#         meanScenario.chance = 1.0
#         meanScenario.pool = []
#         for s in 1:stoc.S
#     		for t in 1:T
#                 meanScenario.data["SL"][t] += stoc.scenarios[s].data["SL"][t]/stoc.S
#     			for b in 1:B
#                     meanScenario.data["SS"][b,t] += stoc.scenarios[s].data["SS"][b,t]/stoc.S
#     			end
#     		end
#         end
#         stoc.scenarios = [meanScenario]
#         stoc.S = 1
#         return stoc
#     end
#
#     function calculate_max_sample(stoc::stocType)
#         meanScenario = scenarioType()
#         meanScenario.ind = 1
#         meanScenario.data = Dict("SL"=>zeros(Float64, T), "SS"=>zeros(Float64, B, T))
#         meanScenario.chance = 1.0
#         meanScenario.pool = []
#         for s in 1:stoc.S
#             for t in 1:T
#                 meanScenario.data["SL"][t] = max(meanScenario.data["SL"][t], stoc.scenarios[s].data["SL"][t])
#                 for b in 1:B
#                     meanScenario.data["SS"][b,t] = max(meanScenario.data["SS"][b,t], stoc.scenarios[s].data["SS"][b,t])
#                 end
#             end
#         end
#         stoc.scenarios = [meanScenario]
#         stoc.S = 1
#         return stoc
#     end
#
#     function calculate_percentil_sample(stoc::stocType, perc::Int)
#         rank = max(floor((perc/100) * stoc.S), 1)
#         percScenario = scenarioType()
#         percScenario.ind = 1
#         percScenario.data = Dict("SL"=>zeros(Float64, T), "SS"=>zeros(Float64, B, T))
#         percScenario.chance = 1.0
#
#         for t in 1:stoc.T
#             sltemp = []
#             sltempnew = []
#             for s in 1:stoc.S
#                 push!(sltempnew, stoc.scenarios[s].data["SL"][t])
#             end
#             percScenario.data["SL"][t] = select(sltempnew, Int(rank))
#
#             for b in 1:B
#                 sstemp_b=[]
#                 sstemp_bnew=[]
#                 for s in 1:stoc.S
#                     push!(sstemp_bnew, stoc.scenarios[s].data["SS"][b,t])
#                 end
#                 percScenario.data["SS"][b,t] = select(sstemp_bnew, Int(rank))
#             end
#         end
#
#         # Refresh the stocType
#         stoc.S = 1
#         stoc.scenarios = [percScenario]
#
#         return stoc
#     end
#
#     function read_stoch_file(filepath::AbstractString)
#
#         stoc = stocType()
#         println("Stochastic file input :: $filepath")
#         if isfile(filepath)
#             stocDict = JSON.parsefile(filepath)
#         elseif isfile(joinpath(config.INPUTPATH, exargs[:PROBLEM], filepath))
#             stocDict = JSON.parsefile(joinpath(config.INPUTPATH, exargs[:PROBLEM], filepath))
#         elseif isfile(joinpath(config.INPUTPATH, filepath))
#             stocDict = JSON.parsefile(joinpath(config.INPUTPATH,filepath))
#         else
#             error("ERROR|stoch.jl|read_stoch_file()|Undetected scenario file. Check your input arguments.")
#         end
#
#         # Validatea file intactness
#         hasSL = false
#         hasSS = false
#         if haskey(stocDict, "SL")
#             hasSL = true
#         else
#             println("No sea level scnearios indicated in this stocFile. Using 0s.")
#         end
#
#         if haskey(stocDict, "SS")
#             hasSS = true
#         else
#             println("No surge scenarios indicate in this stocFile. Using 0s.")
#         end
#
#         if !hasSL && !hasSS
#             error("ERROR|stoch.jl|read_stoch_file()|No scenarios indicated at all")
#         end
#
#         # If There must be at least one scenario
#         if hasSL
#             stoc.S = length(stocDict["SL"])
#             stoc.T = length(stocDict["SL"]["1"])
#             stoc.B = B      #Just use the input value
#             if haskey(stocDict, "S")
#                 @assert stoc.S == stocDict["S"]
#             end
#             if haskey(stocDict, "T")
#                 @assert stoc.T == stocDict["T"]
#             end
#         end
#
#         if hasSS
#             stoc.S = length(stocDict["SS"])
#             stoc.T = length(stocDict["SS"]["1"])
#             stoc.B = length(stocDict["SS"]["1"][1])
#             if haskey(stocDict, "S")
#                 @assert stoc.S == stocDict["S"]
#             end
#             if haskey(stocDict, "T")
#                 @assert stoc.T <= stocDict["T"]
#             end
#             @assert B == stoc.B     # Make sure the stoch file match the current network
#         end
#
#         # Check user defined time steps
#         @assert exargs[:T] == T
#         userT = exargs[:T]
#         if userT > stoc.T
#             error("Trying to run model with high resolution failed. Scenario files is with less.")
#         elseif userT == stoc.T
#             println("User model resolution matchs scenario inputs. No additional action needed.")
#             userSteps = [1:stoc.T;]
#         else
#             println("User model resolution coraser. With fixed horizon... resampling by corasen resolutions.")
#             userSteps = []
#             if stoc.T % userT == 0
#                 userStepSize = Int(stoc.T/userT)
#                 for t in 1:userT
#                     push!(userSteps, t * userStepSize)
#                 end
#                 println("Scenario file input resolution $(stoc.T). Request resolution $(userT). ")
#                 println("Extracting data by steps $(userSteps)")
#                 println("Overriding stoc.T ")
#                 stoc.T = userT
#             else
#                 error("Doesn't support corasion with mismatched data. No approximation on scnearios anymore")
#             end
#         end
#
#         # Start writing scenarios in
#         stoc.scenarios = []
#         stoc.sbdColumns = []
#         for i in 1:stoc.S
#             oneScenario = scenarioType()
#             sl = zeros(Float64,T)
#             ss = zeros(Float64,B,T)
#             j = 1
#             for t in userSteps
#                 if hasSL
#                     sl[j] = stocDict["SL"][string(i)][t]
#                 end
#                 if hasSS
#                     for b in 1:stoc.B
#                         ss[b,j] = stocDict["SS"][string(i)][t][b]
#                     end
#                 end
#                 j += 1
#             end
#             oneScenario.ind = i
#             oneScenario.data = Dict("SL"=>sl, "SS"=>ss)
#             oneScenario.pool = []
#             push!(stoc.scenarios, oneScenario)
#         end
#
#         return stoc
#     end
#
#     # Parsed arg driving
#     splitMode = split(mode, "-")
#     splitLength = length(splitMode)
#     mainMode = ""
#     statMode = ""
#
#     # Detect if serverity is given
#     provideServerity = false
#     mainMode = splitMode[1] #Positional
#     if splitLength > 1
#         if splitMode[2] != "average" && splitMode[2] != "max" && isa(parse(splitMode[2]), Number) != true
#             mainMode = string(splitMode[1],"-",splitMode[2])
#             provideServerity = true
#         end
#     end
#
#     # Detect stat mode
#     provideStat = false
#     if provideServerity == false
#         if splitLength >= 2
#             statMode = splitMode[2]
#             if isa(parse(splitMode[2]), Number) == true
#                 @assert splitMode[3] == "perc"
#             end
#             provideStat = true
#         end
#     else
#         if splitLength >= 3
#             statMode = splitMode[3]
#             if isa(parse(splitMode[3]), Number) == true
#                 @assert splitMode[4] == "perc"
#             end
#             provideStat = true
#         end
#     end
#
#     # Generate main mode
#     if mainMode == "peaceful"
#         samples = generate_zero_sample(samples)
#     elseif mainMode == "random"
#         samples = generate_random_sample(samples, 0.0, 2.0);
#     elseif mainMode == "evolving"
#         samples = generate_evolving_sample(samples, 0.0, 1.0);
#     elseif mainMode == "tuning"
#         samples = generate_tuning_sample(samples);
#     elseif mainMode == "fierce"
#         samples = generate_random_sample(samples, 0.0, 5.0);
#     elseif mainMode == "evolving-highslr"
#         samples = generate_evolving_sample(samples, 3.0, 0.7);
#     elseif mainMode == "evolving-fierce"
#         samples = generate_evolving_sample(samples, 0.0, 1.3);
#     elseif mainMode == "file"
#         if haskey(extras, :filepath)
#             samples = read_stoch_file(extras[:filepath])
#             println("Successfully read stoch file.")
#         else
#             error("Please indicate stoch file path.")
#         end
#     else
#         error("ERROR|stoch.jl|create_samples()|Unkown main mode for sceanrio generation.")
#     end
#
#     if provideStat == true
#         if statMode == "average"
#             samples = calculate_mean_sample(samples)
#         elseif statMode == "max"
#             samples = calculate_max_sample(samples)
#         elseif isa(parse(statMode), Number) == true
#             samples = calculate_percentil_sample(samples, parse(statMode))
#         end
#     end
#
#     # Insert the invisible null scenario
#     dummyScenario = one_null_scenario()
#     push!(samples.scenarios, one_null_scenario())
#
#     # Correct Driver Arguments about the number of samples
#     exargs[:S] = samples.S
#
#     return samples, exargs
# end

function create_null_samples(B::Int, T::Int)

    samples = stocType(1, T, B)
    samples.S = 1
    samples.T = T
    samples.B = B

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

    # println("[STOCH] Total Scenario $(stoc.S)")
    # println("[STOCH] Total Time Step Count $(stoc.T)")
    # println("[STOCH] Data shape => SLR $(size(stoc.scenarios[1].data["SL"]))")
    # println("[STOCH] Data shape => SS $(size(stoc.scenarios[1].data["SS"]))")
    # println("[STOCH] Minimum Bus Elevation => $(minimum(param[:Ele]))")
    # max_slr = 0.0
    # for s in 1:stoc.S
    #     max_slr = max(max_slr, maximum(stoc.scenarios[s].data["SL"]))
    # end
    # println("[STOCH] MAX-SLR => $(max_slr)")
    # println("[STOCH] BUS under MAX-SLR => $(length(param[:Ele][param[:Ele] .<= max_slr]))")
    # for s in 1:stoc.S
    #     println("[STOCH][S=$s] SLR $(stoc.scenarios[s].data["SL"])")
    #     println("[STOCH][S=$s] BUS under SS ")
    #     for t in 1:stoc.T
    #         ss_bus_cnt = 0
    #         for b in 1:param[:B]
    #             if param[:Ele][b] < stoc.scenarios[s].data["SS"][b,t]
    #                 ss_bus_cnt += 1
    #                 println("BUS $(b) surged by $(round(stoc.scenarios[s].data["SS"][b,t]-param[:Ele][b],2))m without protection. (LOAD = $(param[:Pd][b]) | CAP = $(param[:PgUB][b]))")
    #             end
    #         end
    #         println("T$(t)=>$(ss_bus_cnt); ")
    #     end
    #     print("\n")
    # end

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
