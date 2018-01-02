"""
    This is shell function for all parameter generation.
    I can foresee that future parameter generation will mainly relie on reading external file in.
"""
function get_parameters(power::Dict, stoc::stocType, exargs::Dict)

    # Check the drivers in exargs to see how user would like to get parameters
    if isempty(exargs[:PARAMFILE])
        # These internal generation method should only be used for testing
        if exargs[:PROBLEM] == "ieee14"  || exargs[:PROBLEM] == "14"
            param = generate_ieee14_parameters(power,stoc,exargs)
        elseif exargs[:PROBLEM] == "ieee118" || exargs[:PROBLEM] == "ieee118orig" || exargs[:PROBLEM] == "118"
            param = generate_ieee118_parameters(power,stoc,exargs)
        elseif exargs[:PROBLEM] == "3bus"
            param = generate_3bus_parameters(power,stoc,exargs)
        elseif exargs[:PROBLEM] == "4bus"
            param = generate_4bus_parameters(power,stoc,exargs)
        elseif exargs[:PROBLEM] == "nest_case300_ieee.m" || exargs[:PROBLEM] == "ieee300"
            error("ERROR|param.jl|get_parameters()|No implementation for ieee300 ")
        end
    elseif isa(eval(parse(exargs[:PARAMFILE])), Function)
        info("Parameter generator given. Using standarlized generator function.")
        param = eval(parse(exargs[:PARAMFILE]))(power,stoc,exargs)
    else
        param = import_param(power, stoc, exargs)
    end

    return param, exargs
end

"""
    Initialize the parameter package with zeros and power network information
    Initial stochastic information are also determinied here.
    This function should seize the source of generate parameter package.
"""
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

    # # Parse Power systems
    # bus = Dict(Int(bus["index"]) => bus for bus in power["bus"])
    # gen = Dict(Int(gen["index"]) => gen for gen in power["gen"])
    # line = Dict(Int(branch["index"]) => branch for branch in power["branch"])
    line = filter((i, branch) -> branch["br_status"] == 1, line)
    # cost = Dict(Int(gencost["index"]) => gencost for gencost in power["gencost"])

    # Dimension defining parameters
    param[:bus]     = bus
    param[:gen]     = gen
    param[:line]    = line
    # param[:cost]    = cost

    param[:S]       = S   = stoc.S
    param[:T]       = T   = exargs[:T]

    correctEps = exargs[:eps]
    exargs[:eps] = correctEps
    param[:eps]     = exargs[:eps]

    param[:B] = B = bus.count;
    param[:L] = L = line.count;
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
    param[:AngleLimit] = pi;
    param[:AngleShiftLimit] = pi / 180 * exargs[:ANGLESHIFTLambda];
    param[:SHEDLambda] = exargs[:SHEDLambda];
    param[:surplusLoad] = 0.0

    return param
end

"""
    This subroutine reads the parameters of a function and returns all of them as
        a entire dictionary. The key naming convention follows the user's intention.
        Hence, when users want to use these parameters, they can locate them fairly
        easily.
    TODO :: add read parameter subsection when data become available
"""
function generate_3bus_parameters(power::Dict, stoc::stocType, exargs::Dict)

    param = init_parameters(power, stoc, exargs)

    T = param[:T]
    S = param[:S]
    B = param[:B]
    L = param[:L]

    param[:PgUB] = 	[148, 170, 30];
    param[:Hbar] =  [10,  10,  3];
    param[:Pgbar] = [10,  10,  3];
    param[:Pg0]  =	[1,	  1,   0];
    param[:H0]	 = 	[0,   0,   0];
    param[:Ele]  =	[13,  11,  3];
    param[:ProM] =  [5,   5,   4];

    Cg =  [180, 350, 200];
    Ch =  [9999,9999,9999];

    Cg, Ch = lambdaCostTune(Cg, Ch, exargs[:tCostLambda])

    param[:Cg] = [];
    param[:Ch] = [];
    for t in 1:T param[:Cg] = [param[:Cg] ; Cg * 1.0] end
    for t in 1:T param[:Ch] = [param[:Ch] ; Ch * 1.0] end

    for f in exargs[:FEATURES]
        if f == "decreasing-cost"
            info("Decreasing costs over time [1%].")
            param[:Cg] = [];
            for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (0.99)^t] end
            param[:Ch] = [];
            for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (0.99)^t] end
        elseif f == "increasing-cost"
            info("Increasing costs over time [1%].")
            param[:Cg] = [];
            for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (1.01)^t] end
            param[:Ch] = [];
            for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (1.01)^t] end
        end
    end

    for b in 1:B
        param[:Pd][b,1] = param[:bus][b]["pd"]*100;
        param[:Edge][b] = Dict();
        param[:Edge][b]["out"] = [];
        param[:Edge][b]["in"] = [];
    end

    for l in 1:L
        from = param[:line][l]["f_bus"]
        to = param[:line][l]["t_bus"]
        param[:Lcap][from, to] = param[:line][l]["rate_a"] * 100;
        param[:EDGE][from, to] = 1;
        push!(param[:Edge][to]["in"],from)
        push!(param[:Edge][from]["out"],to)
    end

    for t in 2:T
        param[:Pd][:,t] = param[:Pd][:,t-1] * 1.5
    end
    # Pre-determination of the deterministic SLR value
    aslDET = zeros(Bool, B, T, S)
    # Translate stochasticity to deterministic information
    for s in 1:stoc.S
        for i in 1:B
            for t in 1:T
                if stoc.scenarios[s].data["SL"][t] >= param[:Ele][i]
                    aslDET[i,t,s] = false # Offline
                else
                    aslDET[i,t,s] = true # Online
                end
            end
        end
    end

    param = get_aslDet(param, stoc)
    param = get_assDet(param, stoc)
    param = resolve_undersea_load_shift(param, stoc)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    # Sanity check of the parameters before using them
    # Some reshaping work may need to be conducted
    param = check_parameter_intactness(param);

    info("Finish generating parameters for IEEE-14...")
    return param
end

"""
    This subroutine reads the parameters of a function and returns all of them as
        a entire dictionary. The key naming convention follows the user's intention.
        Hence, when users want to use these parameters, they can locate them fairly
        easily.
    TODO :: add read parameter subsection when data become available
"""
function generate_4bus_parameters(power::Dict, stoc::stocType, exargs::Dict)

    param = init_parameters(power, stoc, exargs)

    T = param[:T]
    S = param[:S]
    B = param[:B]
    L = param[:L]

    param[:PgUB] = 	[53,  0,   454, 0]
    param[:Hbar] =  [10,  10,  10,  10]
    param[:Pgbar] = [8,   0,   3,   0]
    param[:Pg0]  =	[1,	  0,   1,   0]
    param[:H0]	 = 	[0,   0,   0,   0]
    param[:Ele]  =	[8,   12,  18,  7]
    param[:ProM] =  [3,   3,   3,   3]

    Cg =  [180, 0,  500, 0]
    Ch =  [30,  30, 80, 50]

    Cg, Ch = lambdaCostTune(Cg, Ch, exargs[:COSTLambda])

    param[:Cg] = [];
    param[:Ch] = [];

    for t in 1:T param[:Cg] = [param[:Cg] ; Cg * 1.0] end
    for t in 1:T param[:Ch] = [param[:Ch] ; Ch * 1.0] end

    for f in exargs[:FEATURES]
        if f == "decreasing-cost"
            info("Decreasing costs over time [1%].")
            param[:Cg] = [];
            for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (0.99)^t] end
            param[:Ch] = [];
            for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (0.99)^t] end
        elseif f == "incresing-cost"
            info("Increasing costs over time [1%].")
            param[:Cg] = [];
            for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (1.01)^t] end
            param[:Ch] = [];
            for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (1.01)^t] end
        end
    end

    for b in 1:B
        param[:Pd][b,1] = param[:bus][b]["pd"] * 100;
        param[:Edge][b] = Dict();
        param[:Edge][b]["out"] = [];
        param[:Edge][b]["in"] = [];
    end

    for l in 1:L
        from = param[:line][l]["f_bus"]
        to = param[:line][l]["t_bus"]
        param[:Lcap][from, to] = param[:line][l]["rate_a"] * 100;
        param[:EDGE][from, to] = 1;
        push!(param[:Edge][to]["in"],from)
        push!(param[:Edge][from]["out"],to)
    end

    for t in 2:T
        param[:Pd][:,t] = param[:Pd][:,t-1] * 1.1
    end

    param = get_aslDet(param, stoc)
    param = get_assDet(param, stoc)
    param = resolve_undersea_load_shift(param, stoc)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    # Sanity check of the parameters before using them
    # Some reshaping work may need to be conducted
    param = check_parameter_intactness(param);

    info("Finish generating parameters for IEEE-14...")
    return param
end

"""
    This subroutine reads the parameters of a function and returns all of them as
        a entire dictionary. The key naming convention follows the user's intention.
        Hence, when users want to use these parameters, they can locate them fairly
        easily.
    TODO :: add read parameter subsection when data become available
    This is the default parameter generator for IEEE-14
"""
function generate_ieee14_parameters(power::Dict, stoc::stocType, exargs::Dict)

    # Initialize a parameter package
    param = init_parameters(power, stoc, exargs)

    T = param[:T]
    S = param[:S]
    B = param[:B]
    L = param[:L]

	# Fill parameter space
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    #                1   2   3   4   5   6   7   8   9   10  11  12  13  14
    param[:PgUB] = 	[208,63, 40, 20, 40, 8,  5,  5,  8,  12, 5,  2,  3,  2];
    param[:Hbar] =  [3,  5,  5,  5,  5,  5,  15, 15, 15, 15, 15, 25, 25, 25];
    param[:Pgbar] = [1,  1,  2,  2,  1,  3,  4,  4,  2,  1,  2,  4,  2,  4];
    param[:Pg0]  =	[1,	 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0];
    param[:H0]	 = 	[0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0];
    param[:Ele]  =	[13, 25, 22, 11, 15, 13, 8,  5,  11, 9,  7,  4,  3,  5];
    param[:ProM] =  [5,  5,  4,  5,  3,  2,  1,  1,  .6, .6, .8, .6, .8, .5];

    # Base costs structure :: stable over time
    Cg =  [180, 350, 200, 130, 160, 150, 70, 40, 100, 80, 50, 20, 15, 20];
    Ch =  [20,  120, 60,  60,  60,  23,  10, 10, 5,  3,  3,  5,  7,  7];

    # Cost tuning for sensitivity analysis
    Cg, Ch = lambdaCostTune(Cg, Ch, exargs[:COSTLambda])

    param[:Cg] = [];
    param[:Ch] = [];
    for t in 1:T param[:Cg] = [param[:Cg] ; Cg * 1.0] end
    for t in 1:T param[:Ch] = [param[:Ch] ; Ch * 1.0] end

    for f in exargs[:FEATURES]
        info("Discount cost over time [$(exargs[:DISCOUNTLambda]*100)%].")
        param[:Cg] = [];
        for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (1 + exargs[:DISCOUNTLambda])^t] end
        param[:Ch] = [];
        for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (1 + exargs[:DISCOUNTLambda])^t] end
    end

    param[:refbus] = 4

    for b in 1:B
        param[:Pd][b,1] = param[:bus][b]["pd"] * 100;
        param[:Edge][b] = Dict();
        param[:Edge][b]["out"] = [];
        param[:Edge][b]["in"] = [];
	end

	for l in 1:L
        from = param[:line][l]["f_bus"]
        to = param[:line][l]["t_bus"]
	    param[:Lcap][from, to] = param[:line][l]["rate_a"] * 100 * exargs[:CONGESTLambda];
	    param[:EDGE][from, to] = 1;
        push!(param[:Edge][to]["in"],from)
        push!(param[:Edge][from]["out"],to)
	end

	for t in 2:T
		param[:Pd][:,t] = param[:Pd][:,t-1] * (1+exargs[:DEMANDLambda])
	end

    # Pre-determination of the deterministic SLR value
    aslDET = zeros(Bool, B, T, S)

    param[:RefBus] = 4

    # Translate stochasticity to deterministic information
    param = get_aslDet(param, stoc)
    param = get_assDet(param, stoc)
    param = resolve_undersea_load_shift(param, stoc)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    # Sanity check of the parameters before using them
    # Some reshaping work may need to be conducted
    param = check_parameter_intactness(param);

    info("Finish generating parameters for IEEE-14...")
    return param
end

"""
    This subroutine reads the parameters of a function and returns all of them as
        a entire dictionary. The key naming convention follows the user's intention.
        Hence, when users want to use these parameters, they can locate them fairly
        easily.
    TODO :: add read parameter subsection when data become available
    This is an available parameter generator for IEEE-118
"""
function generate_ieee118_parameters(power::Dict, stoc::stocType, exargs::Dict)

    param = init_parameters(power, stoc, exargs)

    T = param[:T]
    S = param[:S]
    B = param[:B]
    L = param[:L]
    G = param[:G]

                #  1    2    3    4    5    6    7    8    9    10
    param[:Ele] = [1.0, 1,   1,   6,   4,   1,   4,   8,   4,   1,      #0
                   6,   8,   8,   6,   6,   4,   5,   6,   7,   5,      #1
                   4,   1,   1,   4,   2,   4,   2,   1,   1,   6,      #2
                   2,   3,   6,   3,   1,   2,   2,   3,   1,   2,      #3
                   2,   2,   1,   1,   2,   2,   2,   3,   3,   4,      #4
                   5,   2,   1,   5,   1,   2,   3,   2,   1,   1,      #5
                   2,   3,   1,   2,   2,   3,   2,   3,   2,   2,      #6
                   4,   3,   5,   2,   2,   1,   2,   1,   1,   3,      #7
                   1,   1,   2,   2,   4,   2,   5,   7,   4,   3,      #8
                   5,   5,   5,   1,   3,   2,   1,   2,   5,   2,      #9
                   2,   3,   3,   2,   4,   1,   3,   4,   5,   5,      #10
                   6,   6,   4,   1,   1,   3,   1,   1]                #11

               #  1    2    3    4    5    6    7    8    9    10
    param[:PgUB]=[ 30,  30,  0,   25,  60,  50,  20,  10,  50,  60,     #1
                  10,  287, 40,  0,   0,   16,  20,  24,  10,  0,      #2
                  30,  36,  24,  25,  267, 286, 20,  20,  20,  40,     #3
                  23,  18,  16,  0,   24,  0,   40,  30,  12,  24,     #4
                  12,  12,  12,  20,  4,   8,   10,  0,   232, 16,     #5
                  12,  50,  16,  0,   20,  24,  10,  16,  16,  20,     #6
                  80,  8,   60,  60,  295, 829, 24,  24,  881, 20,     #7
                  20,  18,  0,   12,  12,  16,  10,  24,  18,  445,    #8
                  250, 16,  20,  80,  0,   0,   200, 40,  568, 40,     #9
                  0,   40,  0,   0,   20,  24,  12,  40,  25,  186,    #10
                  30,  30,  98,  24,  150, 16,  40,  120, 120, 200,    #11
                  200, 200, 0,   20,  10,  12,  20,  10]               #12

    #v3           #  1       2        3       4       5       6       7       8       9       10
    #  1       2        3       4       5       6       7       8       9       10
    Cg          = [1126.54,1139.81, 0,      1154.41,1965.47,1139.81,1278.58,1396.21,1636.97,1151.31,
                   1193.01,6501.88, 2139.05,0,      0,      1278.58,2349.28,2305.47,1486.44,0,
                   1636.97,1188.83, 1144.41,1105.06,3616.15,4771.65,1137.14,1126.27,1152.53,1443.42,
                   1371.2, 1120.09, 2023.55,0,      1122.21, 0,     1258.69,1313.54,877.66, 1125.62,
                   1162.81,1125.62, 1122.21,1140.23,1168.12,888.57, 1175.07,1180.00,3729.01,1191.9,
                   1128.43,1145.04, 888.71, 0,      1115.92,1183.75,0,      1132.49,854.59, 1117.41,
                   1392.1, 1115.75, 1192.42,1392.1, 2716.99,5861.28,1251.24,1461.55,5307.96,1144.41,
                   1236.87,1413.73, 0,      1198.49,1198.49,885.74, 1125.99,1129.61,1129.61,5038.12,
                   1728.23,1115.32, 1165.75,1517.38,0,      0,      3685.8, 1737.3, 3595.02,1154.59,
                   0,      1674.64, 0,      0,      1129.84,1125.62,1144.41,1258.69,2349.28,1651.32,
                   1185.34,1156.77, 1813.61,1125.62,3828.97,1111.49,1313.54,3828.97,3953.61,5371.59,
                   6061.14,6061.14, 0,      1139.81,1179.62,1153.85,1179.62,1179.62]
    Cg = Cg * 0.6


               #  1       2       3       4       5       6       7       8       9       10
    Ch        =  [3.77,   3.77,   20,     17.13,  19.89,  3.77,   13.78,  22.84,  15.07,  5.85,
                  17.13,  80.43,  52.48,  20.0,   10.0,   13.78,  24.86,  24.32,  19.99,  15.0,
                  15.07,  4.05,   4.05,   11.42,  28.04,  57.7,   5.71,   2.86,   2.86,   29.83,
                  10.52,  6.98,   20.68,  20.00,  4.05,   10.0,   13.12,  14.92,  2.61,   8.11,
                  8.11,   8.11,   4.05,   4.97,   3.96,   3.96,   7.54,   10.00,  39.76,  10.45,
                  13.06,  7.54,   3.45,   25.0,   3.77,   8.11,   10.00,  5.22,   2.61,   3.45,
                  17.31,  5.94,   8.66,   17.31,  29.19,  66.18,  8.11,   12.16,  45.2,   4.05,
                  16.21,  11.3,   15.0,   6.89,   6.89,   2.61,   5.22,   4.05,   4.05,   51.62,
                  16.48,  2.61,   6.14,   13.12,  15.0,   15.0,   62.44,  26.37,  75.84,  8.57,
                  20.0,   24.86,  10.0,   20.0,   7.83,   8.11,   4.05,   13.12,  24.86,  24.28,
                  9.94,   14.92,  28.16,  8.11,   40.72,  2.61,   14.92,  40.72,  50.9,   62.44,
                  74.93,  74.93,  15.0,   3.77,   3.77,   12.16,  3.77,   3.77]
    Ch = Ch * 6

                  # 1     2     3     4     5     6     7     8     9     10
    param[:ProM] = [0.25, 0.25, 0.25, 1.5,  1,    0.25, 1,    2,    1,    0.25,
                    1.5,  2,    2,    1.5,  1.5,  1,    1.25, 1.5,  1.75, 1.25,
                    1,    0.25, 0.25, 1,    0.5,  1,    0.5,  0.25, 0.25, 1.5,
                    0.5,  0.75, 1.5,  0.75, 0.25, 0.5,  0.5,  0.75, 0.25, 0.5,
                    0.5,  0.5,  0.25, 0.25, 0.5,  0.5,  0.5,  0.75, 0.75, 1,
                    1.25, 0.5,  0.25, 1.25, 0.25, 0.5,  0.75, 0.5,  0.25, 0.25,
                    0.5,  0.75, 0.25, 0.5,  0.5,  0.75, 0.5,  0.75, 0.5,  0.25,
                    1,    0.75, 1.25, 0.5,  0.5,  0.25, 0.5,  0.25, 0.25, 0.75,
                    0.25, 0.25, 0.5,  0.5,  1,    0.5,  1.25, 1.75, 1,    0.75,
                    1.25, 1.25, 1.25, 0.25, 0.75, 0.5,  0.25, 0.5,  1.25, 0.5,
                    0.5,  0.75, 0.75, 0.5,  1,    0.25, 0.75, 1,    1.25, 1.25,
                    1.5,  1.5,  1,    0.25, 0.25, 0.75, 0.25, 0.25 ]

                #   1   2   3   4   5   6   7   8   9   10
    param[:Pgbar]= [3,  2,  0,  3,  1,  2,  2,  3,  1,  3,  #1
                    4,  1,  1,  0,  0,  2,  1,  1,  2,  0,  #2
                    1,  1,  2,  4,  1,  1,  4,  2,  1,  4,  #3
                    1,  8,  1,  0,  4,  0,  2,  2,  6,  2,  #4
                    4,  2,  4,  3,  8,  8,  3,  0,  1,  4,  #5
                    4,  5,  8,  0,  5,  3,  0,  4,  8,  4,  #6
                    2,  8,  3,  2,  1,  1,  1,  1,  1,  2,  #7
                    3,  1,  0,  2,  2,  8,  5,  3,  3,  1,  #8
                    1,  3,  1,  1,  0,  0,  2,  2,  1,  5,  #9
                    0,  2,  0,  0,  8,  2,  2,  2,  1,  1,  #10
                    4,  4,  2,  2,  1,  4,  2,  1,  1,  1,  #11
                    1,  1,  0,  2,  1,  3,  1,  1 ]         #12

                #   1    2    3    4    5    6    7    8    9    10
    param[:Hbar] = [15,  15,  15,  8,   8,   15,  8,   8,   8,   15, #1
                    8,   8,   8,   5,   5,   8,   8,   8,   8,   5,  #2
                    8,   15,  15,  8,   15,  8,   15,  15,  15,  8,  #3
                    15,  8,   8,   5,   15,  5,   15,  8,   15,  15, #4
                    15,  15,  15,  15,  15,  15,  15,  5,   8,   8,  #5
                    8,   15,  15,  5,   15,  15,  5,   15,  15,  15, #6
                    15,  8,   15,  15,  15,  8,   15,  8,   15,  15, #7
                    8,   8,   5,   15,  15,  15,  15,  15,  15,  8,  #8
                    15,  15,  15,  15,  5,   5,   8,   8,   8,   8,  #9
                    5,   8,   5,   5,   8,   15,  15,  15,  8,   15, #10
                    15,  8,   8,   15,  8,   15,  8,   8,   8,   8,  #11
                    8,   8,   5,   15,  15,  8,   15,  15  ]         #12

    param[:Hbar] = param[:Hbar]*10

    Cg, Ch = lambdaCostTune(Cg, Ch, exargs[:COSTLambda])

    param[:Cg] = [];
    param[:Ch] = [];
    for t in 1:T param[:Cg] = [param[:Cg] ; Cg * 1.0] end
    for t in 1:T param[:Ch] = [param[:Ch] ; Ch * 1.0] end

    for f in exargs[:FEATURES]
        info("Discount cost over time [$(exargs[:DISCOUNTLambda]*100)%].")
        param[:Cg] = [];
        for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (1 + exargs[:DISCOUNTLambda])^t] end
        param[:Ch] = [];
        for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (1 + exargs[:DISCOUNTLambda])^t] end
    end

    for g in 1:G
        if param[:gen][g]["pg"] > 0.1
            param[:Pg0][param[:gen][g]["gen_bus"]] = 1.0
        end
    end

    param[:RefBus] = 23
    for b in 1:B
        param[:Ele][b] += rand()  #Lift the buses elevation by a small random portiton to make case more natural
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
			info("Found double coordior ($to, $from)")
            param[:Lcap][from, to] += param[:Lcap][from, to]
		else
			push!(param[:Edge][to]["in"],from)
			push!(param[:Edge][from]["out"],to)
		end
	end

    for t in 2:T
        param[:Pd][:,t] = param[:Pd][:,t-1] * (1+exargs[:DEMANDLambda])
    end
    info("Total Demand Growth $(sum(param[:Pd][:,T])-sum(param[:Pd][:,1])) MW")

    param = get_aslDet(param, stoc)
    param = get_assDet(param, stoc)
    param = resolve_undersea_load_shift(param, stoc)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    param = check_parameter_intactness(param)

    info("Finish generating parameters for IEEE-118...")
    return param
end

"""
    This subroutine reads the parameters of a function and returns all of them as
        a entire dictionary. The key naming convention follows the user's intention.
        Hence, when users want to use these parameters, they can locate them fairly
        easily.
    TODO :: add read parameter subsection when data become available
    This is the default parameter generator for IEEE-118
"""
function paper_param(power::Dict, stoc::stocType, exargs::Dict)

    param = init_parameters(power, stoc, exargs)

    T = param[:T]
    S = param[:S]
    B = param[:B]
    L = param[:L]
    G = param[:G]

    param[:XCoord] = [
            -76.26303523,-76.30822549,-76.31048901,-76.31516578,-76.20887516,-76.25933624,-76.34420679,
            -76.14141991,-76.07308563,-76.01100257,-76.38021394,-76.35656931,-76.43898237,-76.42104946,
            -76.4664604,-76.35629771,-76.39123055,-76.4923841,-76.5341704,-76.59411105,-76.61033868,
            -76.6076774,-76.58130959,-77.09244177,-76.35011927,-76.35338487,-76.17365378,-76.09135172,
            -76.09858182,-76.50251272,-76.26722101,-76.28603685,-76.52047561,-76.44349347,-76.47760819,
            -76.41818152,-76.49812404,-76.62915032,-76.60034086,-76.57867891,-76.58641109,-76.61072064,
            -76.42319383,-76.36450281,-76.60953124,-76.74937004,-76.80276787,-76.70597004,-76.68967515,
            -76.52213105,-76.58488823,-76.54485519,-76.51091048,-76.46601943,-76.44308015,-76.41334089,
            -76.44888633,-76.48486729,-76.46854789,-76.72755167,-76.80467258,-76.79656718,-76.50018542,
            -76.6160793,-76.67525368,-76.7040602,-76.78784714,-76.96435194,-77.1029729,-77.35455688,
            -77.08308332,-77.14968642,-77.0562107,-77.27484889,-77.41016826,-77.39895159,-77.43352044,
            -77.22262804,-77.12102327,-77.25569846,-77.09618294,-77.44460222,-77.43083732,-77.51323306,
            -77.52152246,-77.54785627,-77.58731586,-77.54302105,-77.48835791,-77.50653447,-77.46214226,
            -77.42625038,-77.39827508,-77.39180435,-77.35677223,-77.32931143,-77.31816906,-77.40217444,
            -77.44174161,-77.42498495,-77.48912591,-77.45651493,-77.61916982,-77.47988526,-77.58400787,
            -77.52020721,-77.51743625,-77.56662288,-77.59171051,-77.62245747,-77.60068557,-77.6922786,
            -76.29765723,-76.23757089,-76.19965746,-76.9335021,-76.39670059,-77.32753508]

    param[:YCoord] = [
            36.87113148,36.86649225,36.88604609,36.91971213,36.89522768,36.85315116,36.86243152,36.86693671,
            36.84976986,36.81976193,36.89791814,36.87750193,36.90004917,36.88373572,36.87904801,36.84024022,
            36.82398932,36.87898128,36.94194576,36.85574702,36.79455341,36.76790745,36.72633162,37.02857006,
            36.7411232,36.76143978,36.73888801,36.7839632,36.82221621,36.78386346,36.83277824,36.77267692,
            36.94669698,37.01773155,37.10786354,37.02552363,37.0611691,36.97781932,37.11874016,37.17324602,
            37.20901889,37.22258836,36.98155679,37.08120477,37.17622654,37.28191166,37.37932408,37.24706833,
            37.26460717,37.16099947,37.24453716,37.23051774,37.23590087,37.29374221,37.35256036,37.3021875,
            37.14742403,37.21906836,37.58737356,37.47988297,37.38228193,37.3043588,37.52632356,37.28408538,
            37.1520482,37.15967905,37.23477933,37.22722317,37.46660813,37.19125116,37.21377908,37.18353357,
            37.25966242,37.30375887,37.24818123,37.35538986,37.39389056,37.50929375,37.71922079,37.60574918,
            37.46404283,37.43986343,37.4758825,37.47997834,37.52874562,37.50371869,37.45405318,37.56875573,
            37.54529447,37.57277917,37.55717352,37.53886434,37.51766194,37.57881582,37.60808124,37.54446677,
            37.60176119,37.6494607,37.69041138,37.57730908,37.58758913,37.57025801,37.69224343,37.75543222,
            37.62957482,37.63356903,37.59986311,37.59509995,37.57736499,37.57033404,37.49626608,37.5127317,
            36.78793359,36.81793695,36.78853884,37.04689891,36.86636072,37.37276747]

    param[:Ele] = [
            3.138105869,2.712117672,3.16654563,3.174158573,5.889638424,3.848806143,3.797418356,6.061178684,
            3.703844786,7.868501663,8.644838333,3.703712702,3.905551672,4.167392731,5.934366226,2.451346874,
            4.337753296,5.636176586,6.46437788,22.16767883,20.92450905,15.47028828,15.60319614,33.81655121,
            2.704977274,3.414967775,3.419559956,3.331027508,1.896054864,6.141426086,4.250165462,1.54769063,
            6.355004787,5.744101524,12.83880043,4.524814606,9.872517586,0.676832139,0.256241441,10.29036617,
            9.755583763,23.07926369,8.932640076,2.270322561,6.971250534,16.00845337,31.57906342,19.04461479,
            26.42994499,11.45205784,16.29359245,17.86438751,13.69013691,1.420281053,0.599188864,1.344009638,
            5.721980095,11.01937294,11.07764435,2.524466991,28.01609993,19.96836662,22.9367485,8.168093681,
            6.303514481,10.32908249,8.546386719,33.7030983,37.08203888,42.78454971,2.716136456,36.8248291,
            6.836083412,13.11521626,27.46998024,44.52595901,40.27641678,36.0848732,46.50846481,45.21868515,
            21.02963448,23.43172264,19.3794899,65.19667053,35.68529892,60.50284195,58.61373901,51.86240387,
            60.0736618,72.49128723,61.7067337,7.983781338,43.19146729,57.56072998,46.54642868,50.73565292,
            59.9106636,54.34599686,59.07651138,56.88932419,54.24412918,60.35589218,80.99684906,66.7751236,
            55.27855301,82.74362183,75.71498108,97.76807404,72.74876404,45.13738251,102.416954,113.4431229,
            2.632145643,4.288489819,0.48563766,27.78654671,4.519628048,21.44224358]

                #  1    2    3    4    5    6    7    8    9    10
    param[:PgUB]=[ 30,  30,  0,   25,  60,  50,  20,  10,  50,  60,     #1
                   10,  287, 40,  0,   0,   16,  20,  24,  10,  0,      #2
                   30,  36,  24,  25,  267, 286, 20,  20,  20,  40,     #3
                   23,  18,  16,  0,   24,  0,   40,  30,  12,  24,     #4
                   12,  12,  12,  20,  4,   8,   10,  0,   232, 16,     #5
                   12,  50,  16,  0,   20,  24,  10,  16,  16,  20,     #6
                   80,  8,   60,  60,  295, 829, 24,  24,  881, 20,     #7
                   20,  18,  0,   12,  12,  16,  10,  24,  18,  445,    #8
                   250, 16,  20,  80,  0,   0,   200, 40,  568, 40,     #9
                   0,   40,  0,   0,   20,  24,  12,  40,  25,  186,    #10
                   30,  30,  98,  24,  150, 16,  40,  120, 120, 200,    #11
                   200, 200, 0,   20,  10,  12,  20,  10]               #12

               #  1       2        3       4       5       6       7       8       9       10
   Cg          = [1126.54,1139.81, 0,      1154.41,1965.47,1139.81,1278.58,1396.21,1636.97,1151.31,
                  1193.01,6501.88, 2139.05,0,      0,      1278.58,2349.28,2305.47,1486.44,0,
                  1636.97,1188.83, 1144.41,1105.06,3616.15,4771.65,1137.14,1126.27,1152.53,1443.42,
                  1371.2, 1120.09, 2023.55,0,      1122.21, 0,     1258.69,1313.54,877.66, 1125.62,
                  1162.81,1125.62, 1122.21,1140.23,1168.12,888.57, 1175.07,1180.00,3729.01,1191.9,
                  1128.43,1145.04, 888.71, 0,      1115.92,1183.75,0,      1132.49,854.59, 1117.41,
                  1392.1, 1115.75, 1192.42,1392.1, 2716.99,5861.28,1251.24,1461.55,5307.96,1144.41,
                  1236.87,1413.73, 0,      1198.49,1198.49,885.74, 1125.99,1129.61,1129.61,5038.12,
                  1728.23,1115.32, 1165.75,1517.38,0,      0,      3685.8, 1737.3, 3595.02,1154.59,
                  0,      1674.64, 0,      0,      1129.84,1125.62,1144.41,1258.69,2349.28,1651.32,
                  1185.34,1156.77, 1813.61,1125.62,3828.97,1111.49,1313.54,3828.97,3953.61,5371.59,
                  6061.14,6061.14, 0,      1139.81,1179.62,1153.85,1179.62,1179.62]
    Cg = Cg * 0.8

               #  1       2       3       4       5       6       7       8       9       10
    Ch        =  [3.77,   3.77,   20,     17.13,  19.89,  3.77,   13.78,  22.84,  15.07,  5.85,
                  17.13,  80.43,  52.48,  20.0,   10.0,   13.78,  24.86,  24.32,  19.99,  15.0,
                  15.07,  4.05,   4.05,   11.42,  28.04,  57.7,   5.71,   2.86,   2.86,   29.83,
                  10.52,  6.98,   20.68,  20.00,  4.05,   10.0,   13.12,  14.92,  2.61,   8.11,
                  8.11,   8.11,   4.05,   4.97,   3.96,   3.96,   7.54,   10.00,  39.76,  10.45,
                  13.06,  7.54,   3.45,   25.0,   3.77,   8.11,   10.00,  5.22,   2.61,   3.45,
                  17.31,  5.94,   8.66,   17.31,  29.19,  66.18,  8.11,   12.16,  45.2,   4.05,
                  16.21,  11.3,   15.0,   6.89,   6.89,   2.61,   5.22,   4.05,   4.05,   51.62,
                  16.48,  2.61,   6.14,   13.12,  15.0,   15.0,   62.44,  26.37,  75.84,  8.57,
                  20.0,   24.86,  10.0,   20.0,   7.83,   8.11,   4.05,   13.12,  24.86,  24.28,
                  9.94,   14.92,  28.16,  8.11,   40.72,  2.61,   14.92,  40.72,  50.9,   62.44,
                  74.93,  74.93,  15.0,   3.77,   3.77,   12.16,  3.77,   3.77]
    Ch = Ch * 10 # Fixed adjustment

                  # 1     2     3     4     5     6     7     8     9     10
    param[:ProM] = [0.25, 0.25, 0.25, 1.5,  1,    0.25, 1,    2,    1,    0.25,
                    1.5,  2,    2,    1.5,  1.5,  1,    1.25, 1.5,  1.75, 1.25,
                    1,    0.25, 0.25, 1,    0.5,  1,    0.5,  0.25, 0.25, 1.5,
                    0.5,  0.75, 1.5,  0.75, 0.25, 0.5,  0.5,  0.75, 0.25, 0.5,
                    0.5,  0.5,  0.25, 0.25, 0.5,  0.5,  0.5,  0.75, 0.75, 1,
                    1.25, 0.5,  0.25, 1.25, 0.25, 0.5,  0.75, 0.5,  0.25, 0.25,
                    0.5,  0.75, 0.25, 0.5,  0.5,  0.75, 0.5,  0.75, 0.5,  0.25,
                    1,    0.75, 1.25, 0.5,  0.5,  0.25, 0.5,  0.25, 0.25, 0.75,
                    0.25, 0.25, 0.5,  0.5,  1,    0.5,  1.25, 1.75, 1,    0.75,
                    1.25, 1.25, 1.25, 0.25, 0.75, 0.5,  0.25, 0.5,  1.25, 0.5,
                    0.5,  0.75, 0.75, 0.5,  1,    0.25, 0.75, 1,    1.25, 1.25,
                    1.5,  1.5,  1,    0.25, 0.25, 0.75, 0.25, 0.25 ]

                #   1   2   3   4   5   6   7   8   9   10
    param[:Pgbar]= [3,  2,  0,  3,  1,  2,  2,  3,  1,  3,  #1
                    4,  1,  1,  0,  0,  2,  1,  1,  2,  0,  #2
                    1,  1,  2,  4,  1,  1,  4,  2,  1,  4,  #3
                    1,  8,  1,  0,  4,  0,  2,  2,  6,  2,  #4
                    4,  2,  4,  3,  8,  8,  3,  0,  1,  4,  #5
                    4,  5,  8,  0,  5,  3,  0,  4,  8,  4,  #6
                    2,  8,  3,  2,  1,  1,  1,  1,  1,  2,  #7
                    3,  1,  0,  2,  2,  8,  5,  3,  3,  1,  #8
                    1,  3,  1,  1,  0,  0,  2,  2,  1,  5,  #9
                    0,  2,  0,  0,  8,  2,  2,  2,  1,  1,  #10
                    4,  4,  2,  2,  1,  4,  2,  1,  1,  1,  #11
                    1,  1,  0,  2,  1,  3,  1,  1 ]         #12
    param[:Pgbar] = param[:Pgbar] * 2

                #   1    2    3    4    5    6    7    8    9    10
    param[:Hbar] = [15,  15,  15,  8,   8,   15,  8,   8,   8,   15, #1
                    8,   8,   8,   5,   5,   8,   8,   8,   8,   5,  #2
                    8,   15,  15,  8,   15,  8,   15,  15,  15,  8,  #3
                    15,  8,   8,   5,   15,  5,   15,  8,   15,  15, #4
                    15,  15,  15,  15,  15,  15,  15,  5,   8,   8,  #5
                    8,   15,  15,  5,   15,  15,  5,   15,  15,  15, #6
                    15,  8,   15,  15,  15,  8,   15,  8,   15,  15, #7
                    8,   8,   5,   15,  15,  15,  15,  15,  15,  8,  #8
                    15,  15,  15,  15,  5,   5,   8,   8,   8,   8,  #9
                    5,   8,   5,   5,   8,   15,  15,  15,  8,   15, #10
                    15,  8,   8,   15,  8,   15,  8,   8,   8,   8,  #11
                    8,   8,   5,   15,  15,  8,   15,  15  ]         #12
    param[:Hbar] = param[:Hbar] * 20

    Cg, Ch = lambdaCostTune(Cg, Ch, exargs[:COSTLambda])

    param[:Cg] = [];
    param[:Ch] = [];
    for t in 1:T param[:Cg] = [param[:Cg] ; Cg * 1.0] end
    for t in 1:T param[:Ch] = [param[:Ch] ; Ch * 1.0] end

    for f in exargs[:FEATURES]
        info("Discount cost over time [$(exargs[:DISCOUNTLambda]*100)%].")
        param[:Cg] = [];
        for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (1 + exargs[:DISCOUNTLambda])^t] end
        param[:Ch] = [];
        for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (1 + exargs[:DISCOUNTLambda])^t] end
    end

    for g in 1:G
        if param[:gen][g]["pg"] > 0.1
            param[:Pg0][param[:gen][g]["gen_bus"]] = 1.0
        end
    end

    param[:RefBus] = 23
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
            info("Found double coordior ($to, $from)")
            param[:Lcap][from, to] += param[:Lcap][from, to]
        else
            push!(param[:Edge][to]["in"],from)
            push!(param[:Edge][from]["out"],to)
        end
    end

    for t in 2:T
        param[:Pd][:,t] = param[:Pd][:,t-1] * (1+exargs[:DEMANDLambda])
    end

    info("Total Demand Growth $(sum(param[:Pd][:,T])-sum(param[:Pd][:,1])) MW")

    param = get_aslDet(param, stoc)
    param = get_assDet(param, stoc)
    param = resolve_undersea_load_shift(param, stoc)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    param = check_parameter_intactness(param);

    info("Finish generating parameters for IEEE-118...")
    return param
end

"""
    This subroutine reads the parameters of a function and returns all of them as
        a entire dictionary. The key naming convention follows the user's intention.
        Hence, when users want to use these parameters, they can locate them fairly
        easily.
    TODO :: add read parameter subsection when data become available
    This is the default parameter generator for IEEE-118
"""
function dlT(power::Dict, stoc::stocType, exargs::Dict)

    param = init_parameters(power, stoc, exargs)

    T = param[:T]
    S = param[:S]
    B = param[:B]
    L = param[:L]
    G = param[:G]

    param[:Ele]  = [i for i in exargs[:CSV][:elev_meter_navd88]]
    param[:PgUB] = zeros(B)
    # Cg           = [i for i in exargs[:CSV][:cost_expansion]]
    # Ch           = [i for i in exargs[:CSV][:cost_hardening]]
    Cg           = zeros(B)
    Ch           = zeros(B)
    param[:ProM] = [i for i in exargs[:CSV][:ProM]]
    param[:Pgbar]= zeros(B)
    param[:Hbar] = fill(4, B)

    for f in exargs[:FEATURES]
        info("Discount cost over time [$(exargs[:DISCOUNTLambda]*100)%].")
        param[:Cg] = [];
        for t in 1:T param[:Cg] = [param[:Cg] ; Cg * (1 + exargs[:DISCOUNTLambda])^t] end
        param[:Ch] = [];
        for t in 1:T param[:Ch] = [param[:Ch] ; Ch * (1 + exargs[:DISCOUNTLambda])^t] end
    end

    for g in 1:G
        if param[:gen][g]["pg"] > 0.1
            param[:Pg0][param[:gen][g]["gen_bus"]] = 1.0
            param[:PgUB][param[:gen][g]["gen_bus"]] += param[:gen][g]["pg"]
            param[:Pgbar][param[:gen][g]["gen_bus"]] = 2.0
        end
    end

    param[:RefBus] = 23
    for b in 1:B
        param[:Pd][b,1] = param[:bus][b]["pd"];
        param[:Edge][b] = Dict();
        param[:Edge][b]["out"] = [];
        param[:Edge][b]["in"] = [];
    end

    @show sum(param[:Pd])

    for i in 1:B
        if param[:PgUB][i] > 0.001
            Cg[i] = param[:PgUB][i] * 965 * 1000
        else
            Cg[i] = 0.1
        end
        Ch[i] = (param[:PgUB][i] + param[:Pd][i]) * 500
        (Ch[i] == 0) && (Ch[i] = 99999999)
        @show i, Ch[i], param[:Pd][i], param[:PgUB][i]
    end

    Cg, Ch = lambdaCostTune(Cg, Ch, exargs[:COSTLambda])

    param[:Cg] = [];
    param[:Ch] = [];
    for t in 1:T param[:Cg] = [param[:Cg] ; Cg * 1.0] end
    for t in 1:T param[:Ch] = [param[:Ch] ; Ch * 1.0] end

    for l in 1:L
        from = param[:line][l]["f_bus"]
        to = param[:line][l]["t_bus"]
        param[:Lcap][from, to] = param[:line][l]["rate_a"] * exargs[:CONGESTLambda]
        (param[:Lcap][from, to] == 0.0) && (param[:Lcap][from, to] = 9999)
        param[:EDGE][from, to] = 1
        if (from in param[:Edge][to]["in"])
            # info("Found double coordior ($to, $from)")
            param[:Lcap][from, to] += param[:Lcap][from, to]
        else
            push!(param[:Edge][to]["in"],from)
            push!(param[:Edge][from]["out"],to)
        end
    end

    info("Static demand $(sum(param[:Pd][:,1])) MW")
    param[:Pd][:,1] = param[:Pd][:,1] * (1+exargs[:DEMANDLambda])
    info("T1 demand $(sum(param[:Pd][:,1])) MW")
    for t in 2:T
        param[:Pd][:,t] = param[:Pd][:,t-1] * (1+exargs[:DEMANDLambda])
        info("$(t) demand $(sum(param[:Pd][:,1])) MW")
    end

    param = get_aslDet(param, stoc)
    param = get_assDet(param, stoc)
    param = resolve_undersea_load_shift(param, stoc)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    param = check_parameter_intactness(param);

    info("Finish generating parameters for IEEE-118...")
    return param
end

"""
    This fucntion is designed to param[:Ch]eck the problem's parameter's intactness
    Also, givne different inputs, the parameters are adaptivly reshpaed for
    	for different constraint construction.
    Currently, we only introduce two parameters that need to be adjusted. Later,
    	when we have the parameters structure, the input argument will be utilizing
        that structure.
    Note this this subroutine is subject to param[:Ch]anges when there is a param[:Ch]ange in
        parameter structure. But it handles different dataset.
    This subroutine is going to become useful when reading parameters from file
"""
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

"""
    OUTDATED
    Output currently using parameters to a json file which can be re-utilized for validation and checking.
"""
function write_param(param::Dict)
    paramF = open(string(config.OUTPUTPATH,"/param00.json", "w"))
    write(paramF, JSON.json(paramF))
end

"""
    OUTDATED
    This is the future version of the parameters setup.
    Each parameter should be adjusted out side the code level. It should be an interface that creates such a file.
"""
function import_param(power::Dict, stoc::stocType, exargs::Dict)

    # TODO : This function need to be updated
    paramPath = exargs[:PARAMFILE]
    info("Importing parameters from $paramPath")

    if isfile(paramPath)
        paramS = JSON.parsefile(paramPath)
    elseif isfile(string(config.INPUTPATH,exargs[:PROBLEM],"/",paramPath))
        paramS = JSON.parsefile(string(config.INPUTPATH,exargs[:PROBLEM],"/",paramPath))
    elseif isfile(string(config.INPUTPATH,exargs[:PROBLEM],"/",paramPath,".json"))
        paramS = JSON.parsefile(string(config.INPUTPATH,exargs[:PROBLEM],"/",paramPath,".json"))
    else
        error("ERROR|param.jl|import_param()|Indicated parameter file doesn't exist.")
    end

    param = init_parameters(power, stoc, exargs)

    param[:S] = S = paramS["S"]
    param[:T] = T = paramS["T"]
    param[:eps] = exargs[:eps]

    param[:B]       = B   = param[:bus].count
    @assert B == paramS["bus"].count

    param[:L]       = L   = param[:line].count
    @assert L == paramS["line"].count

    param[:G]       = G   = param[:gen].count;
    @assert G == paramS["gen"].count

    param[:Pd]      = arrarr2mat(paramS["Pd"], Float64, B, T)
    param[:Pgbar]   = paramS["Pgbar"]
    param[:PgUB]    = paramS["PgUB"]
    param[:Pg0]     = paramS["Pg0"]
    param[:H0]      = paramS["H0"]
    param[:Ele]     = paramS["Ele"]
    param[:Hbar]    = paramS["Hbar"]
    param[:ProM]    = paramS["ProM"]

    for b in 1:B
        param[:Edge][b] = Dict();
        param[:Edge][b]["out"] = [];
        param[:Edge][b]["in"] = [];
    end

    for l in 1:L
        from = param[:line][l]["f_bus"]
        to = param[:line][l]["t_bus"]
        param[:Lcap][from, to] = param[:line][l]["rate_a"] * 100;
        param[:EDGE][from, to] = 1;
        push!(param[:Edge][to]["in"],from)
        push!(param[:Edge][from]["out"],to)
    end

    param[:AngleLimit] = paramS["AngleLimit"]
    param[:Cg]      = arrarr2mat(paramS["Cg"],Float64,B,T)
    param[:Ch]      = arrarr2mat(paramS["Ch"],Float64,B,T)
    param[:RefBus]  = paramS["RefBus"]

    param = get_aslDet(param, stoc)
    param = get_assDet(param, stoc)
    param = resolve_undersea_load_shift(param, stoc)

    param = check_parameter_intactness(param);

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
                    @show "SLR Flooding Scenario $(s) :: BUS $(i) || TIME $(t)  (ELE $(param[:Ele][i]), SLR $(stoc.scenarios[s].data["SL"][t]))"
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

"""
    Curretly not in use.
"""
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


"""
    Used for sensitivity test. Balance between expansion costs and harden costs using lambda [0,1].
"""
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

"""
    This function resolves the undersea-load-shift on a higher level.
    Transfering a two dimension (B,T) parameter demand into three dimensions (B, T, S)
    Bad function name
"""
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
