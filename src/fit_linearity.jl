"""
    fit_linearity(pol_order::Int, µ::AbstractVector{<:Union{Real,Measurement{<:Real}}}, peaks::AbstractVector{<:Union{Real,Measurement{<:Real}}}; pull_t::Vector{<:NamedTuple}=fill(NamedTuple(), pol_order+1), v_init::Vector = [], uncertainty::Bool=true )
Fit the calibration lines with polynomial function of pol_order order
    pol_order == 1 -> linear function
    pol_order == 2 -> quadratic function
# Returns
    * `result`: NamedTuple with the following fields
        * `par`: best-fit parameters
        * `gof`: godness of fit
    * `report`: 
"""
function fit_linearity(pol_order::Int, µ_pulser::AbstractVector{<:Measurement{<:T}}, µ_fit::AbstractVector{<:Measurement{<:T}}; e_expression::Union{Symbol, String}="e_ADC", uncertainty::Bool=true) where T<:Real
    @assert length(µ_pulser) == length(μ_fit) "Number of pulser points does not match the number of peak positions"
    @assert pol_order >= 1 "The polynomial order must be greater than 0"

    @debug "Fit calibration curve with $(pol_order)-order polynominal function"
    p_start = append!([0.0, 1.0], fill(0.0, pol_order-1))
    @debug "Initial parameters: $p_start"
    pseudo_prior = _get_fit_calibration_pseudo_prior(pol_order)
    @debug "Pseudo prior: $pseudo_prior"

    # fit calibration curve
    result_fit, report_fit = chi2fit(pol_order, µ_pulser, µ_fit; v_init=p_start, pseudo_prior=pseudo_prior, uncertainty=uncertainty)
    par = result_fit.par

    # built function in string 
    func = join(["$(mvalue(par[i])) * ($(e_expression))^$(i-1)" for i in eachindex(par)], " + ")
    func_err = join(["($(par[i])) * ($(e_expression))^$(i-1)" for i in eachindex(par)], " + ")
    
    result = merge(result_fit, (func = func, func_err = func_err, µ_fit = μ_fit, µ_pulser = µ_pulser))
    report = merge(report_fit, (par = par,))
    return result, report
end
export fit_linearity

function _get_fit_calibration_pseudo_prior(pol_order::Int)
    LegendSpecFits.unshaped(if pol_order == 0
        LegendSpecFits.NamedTupleDist(
            intercept = Normal(0.0, 0.5)
        )
    elseif pol_order == 1
        LegendSpecFits.NamedTupleDist(
            intercept = Normal(0.0, 0.5),
            slope = Normal(1, 0.02)
        )
    elseif pol_order == 2
        LegendSpecFits.NamedTupleDist(
            intercept = Normal(0.0, 0.5),
            slope = Normal(1, 0.02),
            quad = Normal(0.0, (0.005)^2)
        )
    else
        throw(ArgumentError("Only 0, 1, 2 order polynominal calibration is supported"))
    end)
end