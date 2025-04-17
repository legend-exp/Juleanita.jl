"""
"""
function process_pulser_linearity end
export process_pulser_linearity

function process_pulser_linearity(data::LegendData, period::DataPeriod, runs::Vector{<:Int64}, category::Union{Symbol, DataCategory}, channel::ChannelId; reprocess::Bool = true, juleana_logo::Bool = false,  nbins::Int = 50 , rel_cut_fit::T = 0.1) where {T<:Real} 
    e_type = :e_trap;


    # load peakfit results 
    µ_fit = fill(Measurements.measurement(NaN, NaN), length(runs))
    σ_fit = fill(Measurements.measurement(NaN, NaN), length(runs))
    pulser_fit = fill(Measurements.measurement(NaN, NaN), length(runs))
    for (i, run) in enumerate(runs)
        try 
            µ_fit[i]      = data.par[category].rpars.ecal[period, DataRun(run), channel].fit.µ
            σ_fit[i]      = data.par[category].rpars.ecal[period, DataRun(run), channel].fit.σ 
            pulser_fit[i] = data.par[category].rpars.ecal[period, DataRun(run), channel].fit_pulser.µ
        catch 
            @warn "No peak fit results found for run $run"
            continue
        end
    end 
    result, report_lin  = fit_linearity(1, pulser_fit, µ_fit)
end
