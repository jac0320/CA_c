using DataFrames, JSON

function convert(name::AbstractString)

    df_bus = readtable("./$(name)/bus.csv")
    df_gen = readtable("./$(name)/gen.csv")
    df_line = readtable("./$(name)/line.csv")
    df_sub = readtable("./$(name)/sub.csv")
    df_xf = readtable("./$(name)/xf.csv")

    bus = Dict()
    branch = Dict()
    gen = Dict()
    sub_seen = Set()
    busNum2idx = Dict()
    sub2idx = Dict()

    bus_idx = 0
    for i in 1:length(df_bus[:Number])
        if df_bus[:SubNumber][i] in sub_seen
            seen = true
            busNum2idx[df_bus[:Number][i]] = sub2idx[df_bus[:SubNumber][i]]
            @show "mapping BUS $(df_bus[:Number][i]) -> SUB $(df_bus[:SubNumber][i]) -> $(sub2idx[df_bus[:SubNumber][i]])"
        else
            bus_idx += 1
            sub2idx[df_bus[:SubNumber][i]] = bus_idx
            busNum2idx[df_bus[:Number][i]] = bus_idx
            @show "mapping bus_idx $(bus_idx) = BUS $(df_bus[:Number][i]) || SUB $(df_bus[:SubNumber][i])"
            push!(sub_seen, df_bus[:SubNumber][i])
            seen = false
        end
        if seen
            string(df_bus[:LoadMW][i]) == "NA" ? loadMW = 0.0 : loadMW = df_bus[:LoadMW][i]
            bus[string(busNum2idx[df_bus[:Number][i]])]["pd"] += loadMW
        else
            string(df_bus[:LoadMW][i]) == "NA" ? loadMW = 0.0 : loadMW = df_bus[:LoadMW][i]
            bus[string(bus_idx)] = Dict("index"=> bus_idx,
                                "zone"     => 1,
                                "bus_i"    => busNum2idx[df_bus[:Number][i]],
                                "bus_type" => 1,
                                "qd"       => 0.0,
                                "gs"       => 0.0,
                                "bs"       => 0.0,
                                "vmax"     => 0.0,
                                "area"     => 1,
                                "bus_name" => string(df_bus[:Number][i]),
                                "vmin"     => 0.0,
                                "va"       => 0.0,
                                "vm"       => 0.0,
                                "base_kv"  => 0.0,
                                "pd"       =>loadMW)
        end
    end

    gen_idx = 0
    for i in 1:length(df_gen[:BusNum])
        gen_idx += 1
        gen[string(gen_idx)] = Dict("qc1max"     => 0.0,
                                    "model"      => 2,
                                    "startup"    => 0.0,
                                    "qc2max"     => 0.0,
                                    "qg"         => 0.0,
                                    "gen_bus"    => busNum2idx[df_gen[:BusNum][i]],
                                    "ramp_10"    => 0.0,
                                    "mbase"      => 100.0,
                                    "pc2"        => 0.0,
                                    "index"      => gen_idx,
                                    "qmax"       => 0.1,
                                    "pc1"        => 0.0,
                                    "ramp_q"     => 0.0,
                                    "ramp_30"    => 0.0,
                                    "apf"        => 0.0,
                                    "pg"         => df_gen[:MW][i],
                                    "shutdown"   => 0.0,
                                    "ramp_agc"   => 0.0,
                                    "pmax"       => 3.324,
                                    "vg"         => 1.06,
                                    "cost"       => [1.0, 1.0, 0.0],
                                    "gen_status" => 1,
                                    "qmin"       => 0.0,
                                    "qc1min"     => 0.0,
                                    "qc2min"     => 0.0,
                                    "pmin"       => 0.0,
                                    "ncost"      => 3)
    end

    branch_idx = 0
    for i in 1:length(df_line[:BranchDeviceType])
        if busNum2idx[df_line[:BusNumFrom][i]] != busNum2idx[df_line[:BusNumTo][i]]
            df_line[:Circuit][i] == "99" ? lt = 9999 : lt = df_line[:LimitMVAA][i]
            @show df_line[:Circuit][i], lt
            branch_idx += 1
            branch[string(branch_idx)] = Dict("br_r"        => df_line[:R][i],
                                              "rate_a"      => lt,
                                              "shift"       => 0.0,
                                              "br_b"        => df_line[:B][i],
                                              "rate_b"      => 0.0,
                                              "br_x"        => df_line[:X][i],
                                              "rate_c"      => 0.0,
                                              "f_bus"       => busNum2idx[df_line[:BusNumFrom][i]],
                                              "br_status"   => 1,
                                              "t_bus"       => busNum2idx[df_line[:BusNumTo][i]],
                                              "index"       => 1,
                                              "angmin"      => -1.0472,
                                              "angmax"      => 1.0472,
                                              "transformer" => false,
                                              "tap"         => 1.0)
        end
    end

    json_f = open("$(name).json", "w")
    json_d = Dict("name"=>name, "version"=>"1", "baseMVA"=>100, "per_unit" => true, "dcline"=>Dict(),
                  "bus"=>bus,
                  "gen"=>gen,
                  "branch"=>branch)

    JSON.print(json_f, json_d)

    return
end

convert(ARGS[1])
