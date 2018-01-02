# For replicate section
using Glob, DataFrames

filelist = glob("*.out", "/Users/sitew/Dropbox/LANL/Output/climate/paper/REPLICATE/")

outs = Dict(:id=>[], :section=>[], :seed=>[], :scen=>[], :eps=>[], :mode=>[],
            :total=>[], :exp=>[], :hard=>[])

idx = 0
for f in filelist


    if is_windows()
        fname = split(split(f, '\\')[end], '.')[1]
    elseif is_apple() || is_linux()
        fname = split(split(f, '/')[end], '.')[1]
    end

    file = open(f, "r")
    nums = []
    er = false
    for l in readlines(file)
        sl = split(l)
        if ("total" in sl) && ("cost" in sl)
            for i in sl
                k = -1.0
                try
                    k = float(i)
                catch e
                    # do nothing if error
                end
                if k > 0.0
                    push!(nums, k)
                end
            end
        elseif "ERROR:" in sl
            er = true
            break
        end
    end

    if er
        info("Error in $f")
    elseif length(nums) == 3
        idx += 1
        fname_sp = split(fname, '_')
        push!(outs[:id], idx)
        push!(outs[:section], fname_sp[2])
        push!(outs[:seed], fname_sp[3])
        push!(outs[:scen], fname_sp[4])
        push!(outs[:eps], fname_sp[5])
        push!(outs[:mode], fname_sp[6])
        push!(outs[:total], nums[1])
        push!(outs[:exp], nums[2])
        push!(outs[:hard], nums[3])
    else
        info("Unkown case in $f")
    end
end

df = DataFrame()
for header in keys(outs)
    df[header] = outs[header]
end

writetable("outs.csv", df)
