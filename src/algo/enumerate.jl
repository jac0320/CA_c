function enumerator(power::Dict, param::Dict, stoc::stocType, exargs::Dict)

	# The solution enumeration is an very exhaustive way for expreiments
	# Hence it is recommanded to never to over 3
	# Generate SDP pool recursively

	pool = []
	enumerate_generate_pool(0, exargs[:LEVEL], [], stoc.S, pool)
	results = Dict()

	if config.PARALLEL
		scen, obj = pmap((a1,a2,a3,a4,a5)->enumerate_solve_subset(a1,a2,a3,a4,a5),
					[power for s in pool],
					[param for s in pool],
					[stoc for s in pool],
					[exargs for s in pool],
					[s for s in pool])
		for i in 1:length(scen)
			results[scen[i]] = obj[i]
		end
	else
		for s in pool
			scen, obj = enumerate_solve_subset(power, param, stoc, exargs, s)
			results[s] = obj
		end
	end

	return
end

"""
	A small subroutine used by enumerate() that solves a subset of the scnearios
"""
function enumerate_solve_subset(power::Dict, param::Dict, stoc::stocType, exargs::Dict, subset)

	oneJointSubprob = sbd_subprob_formulation(power, param, stoc, subset, exargs)
	solver_config(oneJointSubprob.model, timelimit=config.TIMELIMITIII, mipgap=0.01, showlog=0, focus="optimality", presolve=1, threads=config.WORKERTHREADS)
	status = solve(oneJointSubprob.model, suppress_warnings=true)
	obj = getobjectivevalue(oneJointSubprob.model)

	design = get_design(oneJointSubprob.model)
	totalcost, expandcost, hardencost = get_design_cost(design, param)

	println("[ENUMERATE]Scenario $(subset): The total cost is $(totalcost) = $(expandcost) + $(hardencost)")

	return subset, obj
end

"""
	A small subroutine that generate scenario hybrid pool recursively
"""
function enumerate_generate_pool(l::Int, L::Int, current::Array, S::Int, pool::Array)
	# @show "show up, $l, $current"
	if l == L
		if !(sort(current) in pool)
			push!(pool, sort(current))
		end
		# @show pool
	else
		for i = 1:S
			if !(i in current)
				next = copy([current;i])
			else
				next = copy(current)
			end
			enumerate_generate_pool(l+1, L, next, S, pool)
		end
	end
end
