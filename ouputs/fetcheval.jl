using ArgParse

function parse_commandline_args()
	s = ArgParseSettings()
	@add_arg_table s begin
		"FILENAME"
			help = ".out file"
			arg_type = AbstractString
			default = ""
		"KEY"
			help = "at least one key word provided"
			arg_type = AbstractString
			default = ""
		"LOGIC"
			help = "logic used to fetch line OR/AND"
			arg_type = AbstractString
			default = "OR"
		"--NEXT"
			help = "print the next x lines when found a target line"
			arg_type = Int
			default = 0
		"--K1"
			help = "additional key word 1"
			arg_type = AbstractString
			default = ""
		"--K2"
			help = "addititonal key word 2"
			arg_type = AbstractString
			default = ""
		"--K3"
			help = "addititonal key word 3"
			arg_type = AbstractString
			default = ""
		"--K4"
			help = "addititonal key word 4"
			arg_type = AbstractString
			default = ""
		"--SPLIT"
			help = "splitor"
			arg_type = AbstractString
			default = " "
		"--SHOWLAST"
			help = "only show the target the next x line"
			action = :store_true
		"--VERBOSE"
			help = "slient the utility and only show results if FALSE"
			action = :store_true
	end
	args = parse_args(s)
	return args
end

function main()
	args = parse_commandline_args()
	
	(args["FILENAME"] == "") || (args["KEY"] == "") && error("Need at least a file and a key word")
	!(isfile(args["FILENAME"]) || isfile(string(args["FILENAME"],".out"))) && error("no output file detected")
	!(args["LOGIC"] == "OR") && !(args["LOGIC"] == "AND") && error("illegal input for selecting logic, must be AND/OR")
	
	keywords = [args["KEY"]]
	(args["K1"] != "") && push!(keywords, args["K1"])
	(args["K2"] != "") && push!(keywords, args["K2"])
	(args["K3"] != "") && push!(keywords, args["K3"])
	(args["K4"] != "") && push!(keywords, args["K4"])
	args["VERBOSE"] && println("Detecting keywords $(keywords)")
	args["VERBOSE"] && println("Using logic $(args["LOGIC"])")

	if isfile(args["FILENAME"])
		outf = open(args["FILENAME"], "r")
	else
		outf = open(string(args["FILENAME"], ".out"), "r")
	end
	
	active = 0
	for l in readlines(outf)
		if args["LOGIC"] == "OR"
			checker = [!contains(l, keywords[i]) for i in 1:length(keywords)]
			checker = !prod(checker)
		elseif args["LOGIC"] == "AND"
			checker = [contains(l, keywords[i]) for i in 1:length(keywords)]
			checker = prod(checker)
		end
		checker && (active = 1+args["NEXT"])
		if active > 0
			active -= 1
			args["SHOWLAST"] ? ((active == 0) && println(l)) : println(l)
		end
	end

	close(outf)
	return
end

main()
