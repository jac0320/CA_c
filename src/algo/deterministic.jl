"""
	This is the main function script that defined the climate problem #
	No function/constant definitation here except main()
	Include that libaray that this problem need to use
"""
function deterministic(power::Dict, param::Dict, stoc::stocType, exargs::Dict)

	climate = deterministic_formulation(power, param, stoc, exargs[:MODEL], exargs)
	# warmstart_heuristic(climate, power, param, stoc, exargs)

	# Solve the model deterministically as the entire problem
	info("Sending to solver...")
	numRows = MathProgBase.numlinconstr(climate.model)
	numCols = climate.model.numCols

	info("Solving a problem with $numCols variables and $numRows constraints.")
	info("Targe optimality gap is $(config.OPTGAP)")
	solver_config(climate.model, mipgap=config.OPTGAP, timelimit=config.TIMELIMIT, MIPFocus="optimality", showlog=1, presolve=1, threads=8)

	start_solve = time()
	status = solve(climate.model)
	solve_time = time() - start_solve
	if status == :Infeasible
		info("Original problem infeasible. Getting iis (this might take a while)...")
		print_iis_gurobi(climate.model)
	end

	solution = get_primal_solution(climate)
	design = get_design(climate.model)
	totalcost, expandcost, hardencost = get_design_cost(design, param)
	info(string("The total cost is $(getobjectivevalue(climate.model)) = $expandcost + $hardencost"))

	if !(config.SOLVER == "Cbc")
		lb=getobjbound(climate.model)
		info("Final lower bound is $lb")
		info("Solving wall time is $(solve_time)")
		write_output_files(power, param, stoc, solution, exargs)
	end

	print_design(solution, param)

	flush(STDOUT)

	return climate, solution
end

# This function is the path towards a sketch room where user have the a lot of flexibility
function deterministic_formulation(prob::Dict, param::Dict, stoc::stocType, create_characteristic, exargs::Dict)

	# Basic procedures for generating a formulation
	p = base_formulation(prob, param, stoc, exargs)
	p = create_characteristic(p, exargs)
	p = formulation_features(p, stoc, exargs)
	# p = trim_bounds(p, stoc, exargs)

	return p

end
