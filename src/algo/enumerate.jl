function enumerator(power::Dict, param::Dict, stoc::stocType, driver::Dict)

	pool = [1:param[:S];]
	results = Dict()

	if driver[:PARALLEL]
		pmap((a1,a2,a3,a4)->enu_solve_subset(a1,a2,a3,a4), [param for s in pool], [stoc for s in pool],[driver for s in pool], [[s] for s in pool])
	else
		for s in pool
			enu_solve_subset(param, stoc, driver, [s])
		end
	end

	return
end

"""
	A small subroutine used by enumerate() that solves a subset of the scnearios
"""
function enu_solve_subset(param::Dict, stoc::stocType, driver::Dict, subset)

	sp = build_sp(param, stoc, subset, driver)
	config_solver(sp.model, driver, timelimit=driver[:TIMELIMITIII])
	status = solve(sp.model)

	status == :Infeasible && print_iis_gurobi(sp.model, driver)

	obj = getobjectivevalue(sp.model)
	design = get_design(sp.model)

	totalcost, expandcost, hardencost = get_design_cost(design, param)
	@assert isapprox(obj, totalcost;atol=1e-3)
	info("[ENUMERATE]Scenario $(subset): The total cost is $(totalcost) = $(expandcost) + $(hardencost)")

	return subset, obj
end
