function deterministic(power::Dict, param::Dict, stoc::stocType, driver::Dict)

	climate = basic_problem(power, param, stoc, driver[:MODEL], driver)
	config_solver(climate.model, driver, focus="optimality", threads=16, showlog=driver[:SHOWLOG])

	st = time()
	status = solve(climate.model)
	solve_time = time() - st
	info("Solve time $(solve_time)")

	if status == :Infeasible
		print_iis_gurobi(climate.model) # Infeasible diagnostic
		warn("Model Infeasible.")
		quit()
	end

	solution = get_primal_solution(climate)
	design = get_design(climate.model)
	totalcost, expandcost, hardencost = get_design_cost(design, param)
	lb=climate.model.objBound

	@assert lb <= totalcost
	@assert isapprox(getobjectivevalue(climate.model), totalcost; atol=1e-4)

	info(string("TOTAL $(round(totalcost,4)) || EXP $(round(expandcost,4)) || HARD $(round(hardencost,4)) || LB $(round(lb,4)) || GAP $(driver[:OPTGAP])"))

	write_output_files(power, param, stoc, solution, driver)
	print_design(solution, param)

	return climate, solution
end

# This function is the path towards a sketch room where user have the a lot of flexibility
function basic_problem(prob::Dict, param::Dict, stoc::stocType, complete_formulation, driver::Dict, selection=[])

	isempty(selection) ? selection = [1:param[:S];] : selection = selection
	p = base_formulation(prob, param, stoc, driver)
	complete_formulation(p, param, driver)

	return p
end
