"""
    process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, min_τ::Quantity{T}, max_τ::Quantity{T}, nbins::Int, rel_cut_fit::Real, peak::Symbol, bl_window::ClosedInterval{<:Unitful.Time{<:T}}, tail_window::ClosedInterval{<:Unitful.Time{<:T}}; reprocess::Bool = false) where T <: Real
   
    process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, pz_config::PropDict, bl_window::ClosedInterval{<:Unitful.Time{<:T}}, tail_window::ClosedInterval{<:Unitful.Time{<:T}}; kwargs...) where T <: Real
    
    process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, pz_config::PropDict, dsp_config::DSPConfig; kwargs...)
    
    process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)

Goal: calculate decay time used for pole-zero correction
Inputs:
- `data::LegendData`: LegendData object
- `period::DataPeriod`: data period
- `run::DataRun`: data run
- `category::Union{Symbol, DataCategory}`: data category, e.g. :cal
- `channel::ChannelId`: channel id
- `min_τ::Quantity{T}`: minimum decay time
- `max_τ::Quantity{T}`: maximum decay time
- `nbins::Int`: number of bins for histogram
- `rel_cut_fit::Real`: relative cut for truncated gauss
- `peak::Symbol`: peak to use for decay time calculation. Can also be :all to use all :raw instead of :jlpeaks tier
- `bl_window::ClosedInterval{<:Unitful.Time{<:T}}`: baseline window
- `tail_window::ClosedInterval{<:Unitful.Time{<:T}}`: tail window
Optional:
- `reprocess::Bool`: reprocess the files or not
- `juleana_logo::Bool`: add juleana logo to plots
"""
function process_decaytime end
export process_decaytime

function process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId,
            min_τ::Quantity{T}, max_τ::Quantity{T}, nbins::Int, rel_cut_fit::Real, peak::Symbol,
            bl_window::ClosedInterval{<:Unitful.Time{<:T}}, tail_window::ClosedInterval{<:Unitful.Time{<:T}}; reprocess::Bool = false) where T <: Real
    
    # check if decaytime pars already exist
    pz_file = joinpath(mkpath(data_path(data.par.rpars.pz[period])), "$(string(run)).json")
    if isfile(pz_file) && !reprocess
        @info "Decay time pars (pz) already exist for $category period $period - run $run - channel $channel - you're done!"
        return
    end

    det_ged = _channel2detector(data, channel)
    @debug "Create pars db"
    mkpath(joinpath(data_path(data.par.rpars.pz), string(period)))
    
    if peak == :all 
        data_peak = read_ldata(data, DataTier(:raw), filekeys, channel)
    else
        data_peak  = read_ldata((peak), data, :jlpeaks, category, period, run, channel)
    end 
    wvfs = data_peak.waveform
    decay_times = dsp_decay_times(wvfs, bl_window, tail_window)

    cuts_τ = cut_single_peak(decay_times, min_τ, max_τ,; n_bins=nbins, relative_cut=rel_cut_fit)
    result, report = fit_single_trunc_gauss(decay_times, cuts_τ)
    
    # plot 
    filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
    fig = Figure()
    LegendMakie.lplot!(report, figsize = (600, 430), titlesize = 17, title = get_plottitle(filekey, det_ged, "Decay Time Distribution"), juleana_logo = false, xlabel = "Decay time ($(unit(decay_times[1])))")
    if maximum(abs.(extrema(report.gof.residuals_norm))) > 5.0
        ylim = ceil(Int, maximum(abs.(extrema(report.gof.residuals_norm))))
        ax2 = [ax for ax in fig.content if typeof(ax) <: Makie.Axis][2]
        Makie.ylims!(ax2, -ylim, ylim)
        fig
    end
    savelfig(LegendMakie.lsavefig, fig, data, filekey, channel, :decay_time)

    @info "Found decay time at $(round(u"µs", result.µ, digits=2)) for channel $channel / det $det_ged"
    result_pz = (τ = result.μ, fit = result)
    writelprops(data.par[category].rpars.pz[period], run, PropDict("$channel" => result_pz))
    @info "Saved pars to disk"
    display(fig)
    return fig
end

function process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId,
    pz_config::PropDict, bl_window::ClosedInterval{<:Unitful.Time{<:T}}, tail_window::ClosedInterval{<:Unitful.Time{<:T}}; kwargs...) where T <: Real
   
    process_decaytime(data, period, run, category, channel, pz_config.min_tau, pz_config.max_tau, pz_config.nbins, pz_config.rel_cut_fit, Symbol(pz_config.peak), bl_window, tail_window; kwargs...)
end 

function process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId,
    pz_config::PropDict, dsp_config::DSPConfig; kwargs...)
   
    process_decaytime(data, period, run, category, channel, pz_config.min_tau, pz_config.max_tau, pz_config.nbins, pz_config.rel_cut_fit, Symbol(pz_config.peak), dsp_config.bl_window, dsp_config.tail_window; kwargs...)
end 

function process_decaytime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
    @info "use default values for pz_config and dsp_config from metadata"
   
    # load config: 
    filekeys   = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    pz_config  = dataprod_config(data).dsp(filekeys[1]).pz.default
    dsp_config = DSPConfig(dataprod_config(data).dsp(filekeys[1]).default)

    process_decaytime(data, period, run, category, channel, pz_config.min_tau, pz_config.max_tau, pz_config.nbins, pz_config.rel_cut_fit, Symbol(pz_config.peak), dsp_config.bl_window, dsp_config.tail_window; kwargs...)
end 

