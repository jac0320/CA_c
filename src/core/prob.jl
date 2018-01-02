function read_power(pName::AbstractString)
	power = Dict()
	if pName == "ieee14"
		power = PowerModels.parse_file("./instance/nesta_case14_ieee.m")
		info("Finish reading IEEE-14 power data...")
		return power
	elseif pName == "ieee118"
		power = PowerModels.parse_file("./instance/nesta_case118_ieee.m")
		info("Finish reading local IEEE-118 power data")
		return power
	elseif pName == "ieee118orig"
		power = PowerModels.parse_file("./instance/nesta_case118_ieee_ORIG.m")
		info("Finish reading local IEEE-118 power data - ORIGINAL version")
		return power
	elseif pName == "3bus"
		power = PowerModels.parse_file("./instance/nesta_case3_lmbd.m")
		info("Finish reading local 3-Bus power data")
		return power
	elseif pName == "4bus"
		power = PowerModels.parse_file("./instance/nesta_case4_gs.m")
		info("Finish reading local 3-Bus power data")
		return power
	else
		found = true
		try
			power = PowerModels.parse_file("$(pName)") # User specified problem path, expecting .m file
		catch e
			found = false
			warn("Attempt | Didn't find .m file in path $(pName)")
		end
		found && return power

		found = true
		try
			power = PowerModels.parse_file("./instance/$(pName).m")
		catch e
			found = false
			warn("Attempt 2-m | Didn't find .m file in path ./instance/$(pName).m")
		end
		found && return power

		found = true
		try
			power = PowerModels.parse_json("./instance/$(pName).json")
		catch e
			found = false
			warn("Attempt 2-JSON | Didn't find .json file in path ./instance/$(pName).json")
		end
		found && return power

		found = true
		try
			power = PowerModels.parse_file("./instance/$(pName)")
		catch e
			found = false
			warn("Attemp 3 | Didn't find .m file in path ./instance/$(pName)")
		end
		found && return power

		found = true
		try
			power = PowerModels.parse_file("$(pName).m")
		catch e
			found = false
			warn("Attemp 4 | Didn't find .m file in path in $(pName).m")
		end
		found && return power

		found = true
		try
			power = PowerModels.parse_json("$(pName).json")
		catch e
			found = false
			warn("Attemp 4-JSON | Didn't find .json file in path in $(pName).json")
		end
		found && return power

		found =true
		try
			power = PowerModels.parse_file("$(config.INPUTPATH)$(pName).m")
		catch e
			found = false
			warn("Attempt 5-m | Didn't find .m file in path of $(config.INPUTPATH)$(pName).m")
		end
		found && return power

		found = true
		try
			power = PowerModels.parse_json("$(config.INPUTPATH)$(pName).json")
		catch e
			found = false
			warn("Attempt 5-JSON | Didn't find .json file in path of $(config.INPUTPATH)$(pName).json")
		end
	end
	error("No instance network found.")
end
