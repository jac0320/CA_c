"""
	Input : prob::Dict | param::Dict | stoc::stocType | soln::DesignType | exargs::Dict | ...
	This subroutine perform a basic evaluation on a given design using a user-specific objective
		metric. (Base formulation) and (attach_scenario) are utilized here, which make the formulation
		applied for different characteristics.
	Output : metric::Float64
"""
function eval_formulation(power::Dict, param::Dict, stoc::stocType, scenario, soln, exargs::Dict, objTarget, builtModel=nothing)

	T = exargs[:T]
	B = exargs[:B]

	if objTarget == "feasibility"
		non = nothing
		subprobType = "tight"
	elseif objTarget == "loadshed" || objTarget == "slack" || objTarget == "slackness"
		exargs[:slack] = true
		subprobType = "slackness"
	elseif objTarget == "debug"
		subprobType = "free"
	else
		error("error|formulation.jl|eval_formulation|Unsupport objective type.")
	end

	# Construct a simple subproblem for testing a scenario's user-specific objective
	evalProb = sp_formulation(power, param, stoc, soln, nothing, builtModel)

	# ============================= NEW MODIFIED SECTION ================================ #
	@assert length(scenario) == 1
	scenario = scenario[1]
	evalProb.param[:assDet] = ones(Int, param[:B], param[:T], 1)
	for i in 1:param[:B]
		for t in 1:param[:T]
			expandVal = getupperbound(evalProb.vars[:pg][i,t])
			hardenVal = getupperbound(evalProb.vars[:h][i,t])
			if stoc.scenarios[scenario].data["SS"][i,t] - (param[:ProM][i]*hardenVal + param[:Ele][i]) > 0
				println("Setting Bus $i T=$t to SS=$(stoc.scenarios[scenario].data["SS"][i,t]), PROTECT=$(param[:ProM][i]*hardenVal + param[:Ele][i])")
				evalProb.param[:assDet][i,t,1] = 0
			end
			@assert	haskey(param, :aslDet)
			# setupperbound(evalProb.vars[:a][i,t,1], param[:aslDet][i,t,scenario] * evalProb.param[:assDet][i,t,1])
			# setlowerbound(evalProb.vars[:a][i,t,1], param[:aslDet][i,t,scenario] * evalProb.param[:assDet][i,t,1])
			# setupperbound(evalProb.vars[:ap][i,t,1], param[:aslDet][i,t,scenario] * evalProb.param[:assDet][i,t,1] * expandVal)
			# setlowerbound(evalProb.vars[:ap][i,t,1], param[:aslDet][i,t,scenario] * evalProb.param[:assDet][i,t,1] * expandVal)
			# setupperbound(evalProb.vars[:ah][i,t,1], evalProb.param[:assDet][i,t,1] * hardenVal)
			# setlowerbound(evalProb.vars[:ah][i,t,1], evalProb.param[:assDet][i,t,1] * hardenVal)
			setupperbound(evalProb.vars[:ass][i,t,1], evalProb.param[:assDet][i,t,1])
			setlowerbound(evalProb.vars[:ass][i,t,1], evalProb.param[:assDet][i,t,1])
		end
	end
	# ============================= NEW MODIFIED SECTION ================================ #

	return evalProb
end
