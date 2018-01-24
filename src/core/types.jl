type scenarioType
	ind::Int				# Scenario index
	data::Dict				# More general representation of the stochastic data
	chance::Float64			# Probability of this scenario
	pool::Array 			# Solution pool
	scenarioType(i,data,c) = new(i,data,c,[])
end

type designType
    k           ::Int           # Identifier of a design
	pg			::Any			# Expansion Design Decisions
	h			::Any			# Hardening Design Decisions
	feamap		::Vector 		# Track the feasibility of this design to other scenarios
	cost		::Float64		# Cost of this design
	coverage	::Float64		# Feasibility coverage of this design
	source		::Array			# Scenario pool that results in this design
	lb			::Float64		# Lower bound from thee optimization
	time		::Float64		# Time took to get this design
    active      ::Bool			# Indicator wether this design is active or not
	designType(k::Int, pg::Any, h::Any) = new(k, pg, h, [], Inf, 0.0, [], -Inf, 0.0, true)
	designType(k::Int, pg::Any, h::Any, cost::Float64) = new(k, pg, h, [], cost, 0.0, [], -Inf, 0.0, true)
	designType(k::Int, pg::Any, h::Any, feamap::Vector, cost::Float64, cov::Float64, source::Array, lb::Float64, time::Float64) = new(k, pg, h, feamap, cost, cov, source, lb, time, true)
	designType() = new(0, [], [], [], Inf, 0.0, [], -Inf, 0.0, false)
end

type stocType
	S			::Int					#Total Sample count
	T			::Int					#Total time steps
	B			::Int					#Total bus counts:a useful place to store this number
	scenarios	::Array{scenarioType}	#Each scenario is a scenarioType
	sbdColumns	::Array{designType}					#Stores non-scenario spcific columns
	meta 		::Dict
	stocType(S,T,B) = new(S,T,B,Vector{scenarioType}(S),[],Dict())
end

type solnType
	primal		::Dict
	dual		::Dict
	incumb		::Dict
	solnType() = new()
end

# ============================================================================ #
# JuMP package must be included for this to work
type oneProblem
	vars		::Dict				# Variable Management
	param		::Dict				# Do I need to keep this?
	model		::JuMP.Model		# JuMP optimization model
	status 		::Any				# Optimality Indicator
	oneProblem() = new(Dict(), Dict(),Model(),:None)
end
# ============================================================================ #
