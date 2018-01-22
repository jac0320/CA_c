function deterministic(power::Dict, param::Dict, stoc::stocType, exargs::Dict)

	climate = basic_problem(power, param, stoc, exargs[:MODEL], exargs)

	numRows = MathProgBase.numlinconstr(climate.model)
	numCols = climate.model.numCols
	info("Solving a problem with $numCols variables and $numRows constraints.")
	info("Optimality gap is $(config.OPTGAP)")
	solver_config(climate.model, mipgap=config.OPTGAP, timelimit=config.TIMELIMIT, MIPFocus="optimality", showlog=1, presolve=1, threads=16)

	st = time()
	status = solve(climate.model)
	solve_time = time() - st
	info("Solve time $(solve_time)")
	status == :Infeasible && print_iis_gurobi(climate.model)
	solution = get_primal_solution(climate)
	design = get_design(climate.model)
	totalcost, expandcost, hardencost = get_design_cost(design, param)
	info(string("Cost $(getobjectivevalue(climate.model)) = $expandcost + $hardencost"))

	lb=getobjbound(climate.model)
	info("Final lower bound is $lb")
	info("Solving wall time is $(solve_time)")
	write_output_files(power, param, stoc, solution, exargs)

	print_design(solution, param)

	return climate, solution
end

# This function is the path towards a sketch room where user have the a lot of flexibility
function basic_problem(prob::Dict, param::Dict, stoc::stocType, complete_formulation, exargs::Dict, selection=[])

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	p = base_formulation(prob, param, stoc, exargs)
	complete_formulation(p, param, exargs)
	# p = trim_bounds(p, stoc, exargs)

	return p
end
