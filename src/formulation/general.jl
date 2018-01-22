function base_formulation(prob::Dict, param::Dict, stoc::stocType, exargs::Dict, selection=[]; kwargs...)

	base = oneProblem()
	solver_config(base.model)
	param = check_parameter_intactness(param)

	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	# Variables :: ========================================================== #
	post_adaptation_vars(base, param)
	post_logical_vars(base, param, selection)

	# Constraints ============================================================ #
	post_incremental_cons(base, param)
	post_logical_cons(base, param, stoc, selection)
	post_risk_cons(base, param, exargs, selection)

	# Objective
	post_adaptation_obj(base, param)

	return base
end

function cb_model(prob::oneProblem, param::Dict, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)
	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	post_cb_vars(prob, param, selection)
	post_cb_cons(prob, param, selection, subprobType)

	return prob
end

function cnf_model(prob::oneProblem, param::Dict, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)
	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	haskey(options, :logical) ? logical = options[:logical] : logical = true
	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	post_cnf_vars(prob, param, selection, subprobType)
	post_cnf_cons(prob, param, selection, subprobType)

	return prob
end

function dcpf_model(prob::oneProblem, param::Dict, exargs::Dict, selection=[]; kwargs...)

	options = Dict(kwargs)

	haskey(options, :subprobType) ? subprobType = options[:subprobType] : subprobType = "free"
	isempty(selection) ? selection = [1:param[:S];] : selection = selection

	B = param[:B]

	# [TODO] this really should be happening here
	E = length(param[:line])
	lineX = zeros(Float64, B, B)
	for l in 1:E
		lineX[param[:line][l]["f_bus"], param[:line][l]["t_bus"]] = param[:line][l]["br_x"]
	end
	param[:lineX] = lineX / 100

	post_dcpf_vars(prob, param, selection, subprobType)
	post_dcpf_cons(prob, param, selection, subprobType)

	return prob
end

function mccormick(m,xy,x,y,xˡ,xᵘ,yˡ,yᵘ)
    @constraint(m, xy >= xˡ*y + yˡ*x - xˡ*yˡ)
    @constraint(m, xy >= xᵘ*y + yᵘ*x - xᵘ*yᵘ)
    @constraint(m, xy <= xˡ*y + yᵘ*x - xˡ*yᵘ)
    @constraint(m, xy <= xᵘ*y + yˡ*x - xᵘ*yˡ)
    return m
end
