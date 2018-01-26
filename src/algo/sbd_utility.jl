"""
	Used to be inside function lineup_cg_heuristic
"""
function lineup_join_two_scenario(power::Dict, param::Dict, stoc::stocType, exargs::Dict, s::Int, solved::Dict, pickScenarioPool,subprob_formulation::Function)

	# Available Public Variable: all upper-level input arguments
	designPool = Array{designType}(0)

	if config.WARMSTART
		allSubprobs = Array{oneProblem}(stoc.S)
		# Prepare a vector of problems for repetitvely checking
		for ss = 1:stoc.S
			allSubprobs[ss] = oneProblem()
			allSubprobs[ss] = sbd_base_formulation(power, param, stoc)
			allSubprobs[ss] = attach_scenario(allSubprobs[ss], stoc, [ss], exargs[:MODEL], 0.0, exargs, subprobType="tight")
		end
	end

	for subs in pickScenarioPool
		if subs != s && stoc.sbdColumns[s].feamap[subs] == 0
			# Mark the problem to be solved
			testPool = [s, subs]
			if haskey(solved, testPool) || haskey(solved, flipdim(testPool,1))
				# info("($(myid()))[Solved] ~ $testPool")
			else
				oneJointSubprob = subprob_formulation(power, param, stoc, testPool, exargs)
				solver_config(oneJointSubprob.model, license=1, timelimit=config.TIMELIMITIII, mipgap=config.OPTGAP, showlog=config.SHOWLOG, focus="optimality", presolve=1, threads=config.WORKERTHREADS)

				status = solve(oneJointSubprob.model, suppress_warnings=true)
				if status == :Infeasible
					print_iis_gurobi(oneJointSubprob.model)
					error("ERROR|alg_sbd_jl|sbd_heuristic|Subprob infeasible, check formulation")
				end

				if config.WARMSTART
					collectTime = @elapsed oneDownwardDesign = umn(power, param, stoc, exargs,
														solved,
														oneJointSubprob,
														testPool,
														builtModel = allSubprobs)
				else
					collectTime = @elapsed oneDownwardDesign = umn(power, param, stoc, exargs,
														solved,
														oneJointSubprob,
														testPool)

				end
				push!(designPool, oneDownwardDesign)
				d = oneDownwardDesign
				info("($(myid()))[Spinning-S $s] $(d.source) cost = <$(d.cost)($(round.(d.coverage*100,2))\%)> [LB >> $(round.(d.lb,2))][Time $(round.(d.time,2))s][Collector $(collectTime)];")
				flush(STDOUT)
			end
		end
	end
	return designPool
end

function lineup_generate_pool(stoc::stocType, setA, setB, solved::Dict)

	probPool = []
	for sa in setA
		for sb in setB
			if sb != sa && stoc.sbdColumns[sa].feamap[sb] == 0 && !haskey(solved, [sa;sb]) && !haskey(solved, [sb;sa])
				push!(probPool, [sa; sb])
			end
		end
	end

	return probPool
end


function partial_generate_pool(stoc::stocType, idxA, setA, solved::Dict)

	probPool = []
	for sa in setA
		if idxA != sa && stoc.sbdColumns[sa].feamap[idxA] == 0 && !haskey(solved, [sa;idxA]) && !haskey(solved, [idxA;sa])
			push!(probPool, [sa; idxA])
		end
	end

	return probPool
end

function lineup_solve_pair(power::Dict, param::Dict, stoc::stocType, exargs::Dict, pair, subprob_formulation)

	onePairSubprob = subprob_formulation(power, param, stoc, pair, exargs)
	solver_config(onePairSubprob.model, license=1, timelimit=config.TIMELIMITIII, mipgap=config.OPTGAP, showlog=config.SHOWLOG, focus="optimality", presolve=1, threads=config.WORKERTHREADS)
	status = solve(onePairSubprob.model, suppress_warnings=true)
	if status == :Infeasible
		print_iis_gurobi(onePairSubprob.model)
		error("ERROR|alg_sbd_jl|sbd_heuristic|Subprob infeasible, check formulation")
	end

	onePairDesign = sbd_collect_column(power, param, stoc, exargs,
										Dict(),
										onePairSubprob,
										pair,
										skipFea = 1)

	return onePairDesign
end

function lineup_regroup_design(designPool)

	designCnt = length(designPool)
	if length(config.WORKERS) > 1
		workerCnt = length(config.WORKERS) * Int(config.PARALLEL) + (1 - Int(config.PARALLEL))
	else
		workerCnt = config.WORKERS[1] * Int(config.PARALLEL) + (1 - Int(config.PARALLEL))
	end

	perGroupCnt = floor(designCnt / workerCnt)
	groupDesigns = []

	if designCnt <= workerCnt
		# Getting ready for testing
		for i in 1:designCnt
			oneGroup = []
			push!(oneGroup, designPool[i])
			push!(groupDesigns, oneGroup)
		end
	else
		k = 0
		for i in 1:workerCnt
			if i != workerCnt
				oneGroup = []
				for j in 1:perGroupCnt
					k += 1
					push!(oneGroup, designPool[k])
				end
			else
				lastGroupCnt = designCnt - k
				oneGroup = []
				for j in 1:lastGroupCnt
					k += 1
					push!(oneGroup, designPool[k])
				end
			end
			push!(groupDesigns, oneGroup)
		end
	end

	return groupDesigns
end

function lineup_groupdesign_fea(power::Dict, param::Dict, stoc::stocType, exargs::Dict, designGroup)

	allSubprobs = Array{oneProblem}(stoc.S)

	# Prepare a vector of problems for repetitvely checking
	for ss = 1:stoc.S
		allSubprobs[ss] = oneProblem()
		allSubprobs[ss] = sbd_base_formulation(power, param, stoc)
		allSubprobs[ss] = attach_scenario(allSubprobs[ss], stoc, [ss], exargs[:MODEL], 0.0, exargs, subprobType="tight")
	end

	for d in designGroup
		#info("[WORKER$(myid())]Start checking feasibility $(now())")
		non,d.coverage,ifP,fP=check_feasible(power,param,stoc,d,exargs,[],precheck=false,builtModel=allSubprobs)
		for subs in fP
			d.feamap[subs] = 1
		end
	end

	return designGroup
end

function sbd_checkLB(design::designType, stoc::stocType, UB::Float64, LB::Float64, exargs::Dict; kwargs...)

	S = stoc.S

	if design.coverage >= (1-exargs[:eps])
		print(">! [", round.(design.cost,4) ,"]",design.source," !<")
		kicker = 0
		if design.cost < UB
			kicker += 1
			for subs in [1:S;]
				if design.cost < stoc.sbdColumns[subs].cost
					kicker += 1
					if !(subs in kickSet)
						push!(kickSet, subs)
					end
				end
			end

			if kicker >= S * exargs[:eps]
				# Calculate gap to fixedLB
				gaptofixedLB = round.((design.cost - LB) / design.cost * 100, 2)
				# This mean we find a new incumbent, not necessarily optimal answer
				print(">![", round.(gaptofixedLB,4), "<!")
				if gaptofixedLB - 0.0 <= config.TOLERANCE
					info("Locating optimal solution... But let's not not stop for now...")
					return true
				end
			end
		end
	end

	return false
end

"""
	Needs update
"""
function sbd_analysis(power::Dict, param::Dict, stoc::stocType, exargs::Dict, onescenario_formulation::Function; kwargs...)
	lockInference = lockbus_inference(power,param,stoc,exargs,onescenario_formulation)
	return 0
end

"""
	Needs update
"""
function lockbus_inference(power::Dict,
							param::Dict,
							stoc::stocType,
							exargs::Dict,
							onescenario_formulation::Function;
							kwargs...)

	B = param[:B]
	T = param[:T]
	S = stoc.S

	# Feasibility Checking
	subprobType = "tight"

	criticalBusCuts = Dict()
	criticalBusCuts[:expand] = Dict()
	criticalBusCuts[:harden] = Dict()

	# Perform expansion analyses on each bus at all scenarios
	for s in 1:S
		for b in 1:B
			oneScenProb = onescenario_formulation(power,param,stoc,[s],exargs,subprobType=subprobType)
			# Fix the expansion upper bound
			for t in 1:T
				setupperbound(oneScenProb.vars[:pg][b,t], param[:Pg0][b])
			end
			status = solve(oneScenProb.model, suppress_warnings=true)
			if status == :Optimal
				# Alternative solution found, able to push further from here (Alternative branch)
			elseif status == :Infeasible
				if haskey(criticalBusCuts[:expand], s)
					push!(criticalBusCuts[:expand][s], [b])
				else
					criticalBusCuts[:expand][s] = []
					push!(criticalBusCuts[:expand][s], [b])
				end
			else
				error("ERROR|sbd.jl|lockbus_inference()|Unkown model status.")
			end
			println("Scenario $s :: Bus $b lock found $status alternative")
		end
	end


	# Similar analyses can be conducted here with hardening
	for s in 1:S
		for b in 1:B
			oneScenProb = onescenario_formulation(power,param,stoc,[s],exargs,subprobType=subprobType)
			# Fix the expansion upper bound
			for t in 1:T
				setupperbound(oneScenProb.vars[:h][b,t], param[:H0][b])
			end
			status = solve(oneScenProb.model, suppress_warnings=true)
			if status == :Optimal
				# Alternative solution found, able to push further from here (Alternative branch)
			elseif status == :Infeasible
				if haskey(criticalBusCuts[:harden], s)
					push!(criticalBusCuts[:harden][s], [b])
				else
					criticalBusCuts[:harden][s] = []
					push!(criticalBusCuts[:harden][s], [b])
				end
			else
				error("ERROR|sbd.jl|lockbus_inference()|Unkown model status.")
			end
			println("Scenario $s :: Bus $b lock found $status alternative")
		end
	end

	return criticalBusCuts
end

# This subroutine collects a design into a scenario's design pool without introducing duplicates in the scenario design pool
function sbd_store_design(designPool::Array, design::designType; kwargs...)

	options = Dict(kwargs)

	if haskey(options, :incumbent)
		incumbent = options[:incumbent]
	else
		incumbent = Inf
	end

	poolLength = length(designPool)
    same = false
	B, T = size(design.pg)

	if design.cost <= incumbent

		skip = false
		for i in 1:poolLength
			if abs(design.cost - incumbent) <= 0.1
				skip = true # Skip the expensive columns
				break
			end
		end

		if !skip
			for i in 1:poolLength
				if (design.pg == designPool[i].pg) && (design.h == designPool[i].h)
					same = true
					break
				end
			end
		end

		# Only insert if no replication has detected
		if !same && !skip
			push!(designPool, design)
		end
	end


    return designPool
end

function sbd_collect_column(power::Dict, param::Dict, stoc::stocType, exargs::Dict, solved::Dict, subprob::oneProblem, sourcePool::Array; kwargs...)

	options = Dict(kwargs)

	haskey(options, :skipFea) ? checkFea = false : checkFea = true	# By default, this is alwasy true

	# Fetch design from model
	design = get_design(subprob.model)
	design.active = true
	design.source = sourcePool
	design.cost = getobjectivevalue(subprob.model)
	design.lb = solver_lower_bound(subprob.model)
	design.time = getsolvetime(subprob.model)
	design.feamap = zeros(Int, exargs[:S])

	if checkFea
		if haskey(options, :builtModel)
			non, design.coverage, ifP, fP = check_feasible(power, param, stoc, design, exargs, [], precheck=false, builtModel = options[:builtModel], licnese="distributed")
		else
			non, design.coverage, ifP, fP = check_feasible(power, param, stoc, design, exargs, [], precheck=false)
		end

		for subs in fP
			design.feamap[subs] = 1
		end
	else
		design.coverage = -1.0
	end

	solved[sourcePool] = design
	solved[flipdim(sourcePool, 1)] = design

	return design
end


"""
	Takes a subset of scenarios and solve using solver or sbd_norisk. Perform a simple optimality check afterward.
	Indexing matters here.
"""
function sbd_ubPool_solve(power::Dict, param::Dict, stoc::stocType, exargs::Dict,
							fixedLB::Float64, ubPool::Array, ubISOCost::Array, dsp_formulation::Function, timeleft=config.TIMELIMITII)

	S = exargs[:S]

	# Generate DSP problem formulation
	if !config.USESBDNORISK
		ubPoolProb = dsp_formulation(power, param, stoc, ubPool, exargs)
		warmstart_heuristic(ubPoolProb, power, param, stoc, exargs, selection=ubPool)
		info("[SBD] UB Pool todo-> $ubPool")
		solver_config(ubPoolProb.model,
						timelimit=timeleft,
						mipgap=config.OPTGAP,
						showlog=config.SHOWLOG,
						focus="optimality",
						presolve=1,
						threads=16)
		status = solve(ubPoolProb.model)
		ubPoolDesign = sbd_collect_column(power, param, stoc, exargs, Dict(), ubPoolProb, ubPool)
		stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, ubPoolDesign, incumbent = Inf)
	else
		ubPoolProb, ubPoolDesign = sbd_norisk(power, param, stoc, exargs, dsp_formulation, ubPool, ubISOCost)
		ubPoolDesign = sbd_collect_column(power, param, stoc, exargs, Dict(), ubPoolProb, ubPool)
		stoc.sbdColumns = sbd_store_design(stoc.sbdColumns, ubPoolDesign, incumbent = Inf)
	end

	info("[SBD] UB pool is solved with objective >> $(ubPoolDesign.cost) [$(ubPoolDesign.coverage)] in [ $(getsolvetime(ubPoolProb.model))s]")

	if ubPoolDesign.cost - fixedLB > config.TOLERANCE
			info("[SBD] UB pool provide a feasible solution.")
	else ubPoolDesign.cost - fixedLB <= config.TOLERANCE
			info("[SBD] UB pool sufficient for optimality.")
		return true, ubPoolProb
	end

	return false, ubPoolProb
end

"""
	This small routine performs an iteration check in all stored columns to eliminate scenarios on second level.
"""
function columns_manager(stoc::stocType, kickScen, incumbent::Float64, exargs::Dict; kwargs...)

	allScen = Set()
	for s in 1:stoc.S
		push!(allScen, s)
	end

	if isa(kickScen, Array)
		kickScenSet = Set(kickScen)
	elseif isa(kickScen, Set)
		kickScenSet = kickScen
	end
	activeScen = setdiff(allScen, kickScenSet)

	eliminateScen = Set()
	scenEliminater = Dict(i=>Set() for i=1:stoc.S)

	S = stoc.S - length(kickScen)
	requireS = ceil(stoc.S * (1-exargs[:eps]))

	# First Summarize All Generated Columns
	for i in 1:length(stoc.sbdColumns)
		if incumbent + config.TOLERANCE < stoc.sbdColumns[i].lb && (length(S) < stoc.S) && length(stoc.sbdColumns[i].source) <= 2
			if stoc.sbdColumns[i].active
				println("[COLUMN MANAGER] Found eliminating column [scenario set $(stoc.sbdColumns[i].source)][Cost $(stoc.sbdColumns[i].cost)]\n")
			end
			stoc.sbdColumns[i].active = false
			# Inference phase
			for s in stoc.sbdColumns[i].source
				avoidSet = setdiff(Set(stoc.sbdColumns[i].source),s)
				if !isempty(avoidSet)
					union(scenEliminater[s], unique(avoidSet))
				end
			end
		end
	end

	for i in activeScen
		combineAvoid = collect(unique(scenEliminater[i]))
		if length(combineAvoid) > stoc.S - requireS
			push!(eliminateScen, i)
		end
	end

	if !isempty(eliminateScen)
		info("[COLUMN MANAGER] Helped eliminate $(length(eliminateScen))")
	end

	return collect(eliminateScen)

end

"""
	This small routine checks all designs and deactivate the ones that are useless
"""
function columns_cleaner(stoc::stocType, incumbent::Float64; kwargs...)

	for i in 1:length(stoc.sbdColumns)
		if stoc.sbdColumns[i].lb >= incumbent + config.TOLERANCE && stoc.sbdColumns[i].active
			info("Deactivating column $(i) with cost $(round.(stoc.sbdColumns[i].cost,2)). [$(round.(incumbent,2))]\n")
			stoc.sbdColumns[i].active = false
		end
	end

	return 0
end

"""
	Print summary information line by line about all collected columns
"""
function columns_printer(stoc::stocType)
	for i in 1:length(stoc.sbdColumns)
		column = string("[COLUMNS] ($i) [COST $(round.(stoc.sbdColumns[i].cost,2))][ACTIVE $(stoc.sbdColumns[i].active)][Cover $(stoc.sbdColumns[i].coverage)]\n")
	end
	return
end

function sbd_solve_master(power::Dict, param::Dict, stoc::stocType, exargs::Dict; kwargs...)

	options = Dict(kwargs)

	if haskey(options, :prevSolution)
		prevSolution = options[:prevSolution]
		@assert isa(prevSolution, solnType)
	else
		prevSolution = []
	end

	if haskey(options, :incumbent)
		master = sbd_master_formulation(power, param, stoc, exargs, exargs[:eps], [], prevSolution, incumbent=options[:incumbent])
	else
		master = sbd_master_formulation(power, param, stoc, exargs, exargs[:eps], [], prevSolution)
	end

	solver_config(master.model,timelimit=config.TIMELIMITII, mipgap=config.OPTGAP, showlog=1, focus="optimality", presolve=1, threads=16)
	status = solve(master.model)

	if status == :Infeasible
		print_iis_gurobi(master.model)
		error("ERROR|alg_sbd.jl|sbd_heuristic()|master proble infeasible, check formulation");
	else
		masterDesign = get_design(master.model)
		masterDesign.lb = solver_lower_bound(master.model)
		masterDesign.time = getsolvetime(master.model)
		masterDesign.cost = getobjectivevalue(master.model)
		non, masterDesign.coverage, non, feaPool = check_feasible(power, param, stoc, masterDesign, exargs, [], precheck=true)
		info("[MASTER] Master Objective = $(masterDesign.cost)[LB=$(masterDesign.lb)][$(masterDesign.time)s][", round.(masterDesign.coverage*100,2), "]")
		info("[MASTER] Master solution information summary...")
		for y in 1:length(master.vars[:Y])
			YValue = getvalue(master.vars[:Y])
			if YValue[y] == 1
				info("[MASTER] Column $y selected : [COST $(stoc.sbdColumns[y].cost)][SOURCE $(stoc.sbdColumns[y].source)][COVER $(stoc.sbdColumns[y].coverage)][TIME $(stoc.sbdColumns[y].time)][LB $(stoc.sbdColumns[y].lb)]")
			end
		end
		masterSolution = get_primal_solution(master)
	end

	# Find out what scenarios did the master problem picked
	pickScenarioPool, neglectedScenarioPool = get_master_selection(master, exargs)
	masterDesign.source = pickScenarioPool

	return masterDesign, pickScenarioPool, neglectedScenarioPool, masterSolution
end
