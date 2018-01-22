module ADCCTypes

	using JuMP
	using ArgParse
	using PowerModels

	export designType, scenarioType, stocType, oneProblem, solnType

	type scenarioType
		ind::Int				# Scenario index
		data::Dict				# More general representation of the stochastic data
		chance::Float64			# Probability of this scenario
		pool::Array 			# Solution pool
		scenarioType(i,data,c) = new(i,data,c,[])
	end

	type stocType
		S			::Int					#Total Sample count
		T			::Int					#Total time steps
		B			::Int					#Total bus counts:a useful place to store this number
		scenarios	::Array{scenarioType}	#Each scenario is a scenarioType
		sbdColumns	::Array					#Stores non-scenario spcific columns
		meta 		::Dict
		stocType(S,T,B) = new(S,T,B,Vector{scenarioType}(S),[],Dict())
	end

	type designType
        k           ::Int           # Identifier of a design
		pg			::Any			# Expansion Design Decisions
		h			::Any			# Hardening Design Decisions
		feamap		::Array 		# Track the feasibility of this design to other scenarios
		cost		::Float64		# Cost of this design
		coverage	::Float64		# Feasibility coverage of this design
		source		::Array			# Scenario pool that results in this design
		lb			::Float64		# Lower bound from thee optimization
		time		::Float64		# Time took to get this design
        active      ::Bool			# Indicator wether this design is active or not
		designType() = new()
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

	type cellType
		name		::AbstractString
		master		::oneProblem
		subprob 	::oneProblem
		stage		::Int
		cellType() = new()
	end

end # Module
