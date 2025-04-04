"""
    progress_peakfits(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; reprocess::Bool = true, e_types::Vector{<:Symbol} = [:e_trap, :e_cusp, :e_zac], juleana_logo::Bool = false)
    Process the peak fits for the benchtest data. 
"""
function process_peakfits end
export process_peakfits
function process_peakfits(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId;
    reprocess::Bool = true, juleana_logo::Bool = false,  nbins::Int = 50 , rel_cut_fit::T = 0.1) where {T<:Real}
    e_type = :e_trap; 

    if !reprocess && haskey(data.par[category].rpars.ecal[period, run, channel], e_type)
        @info "Energy calibration already exists for all $(e_types)  -> you're done!"
        return
    end

    # load plot folder and detector
    filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
    det = _channel2detector(data, channel)
    plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekey, :peakfits) * "/"
    
    # load energy estimates from dsp pars 
    dsp_pars = read_ldata(data, :jldsp, category, period, run, channel);
  
    # waveforms 
    e_uncal = filter!(isfinite, dsp_pars[e_type])
    qmin, qmax = 0.0, 1.0
    emin = quantile(e_uncal, qmin)
    emax = quantile(e_uncal, qmax)
    cuts_τ = cut_single_peak(e_uncal, emin, emax,; n_bins=nbins, relative_cut=rel_cut_fit)
    result, report = fit_single_trunc_gauss(e_uncal, cuts_τ)
    @debug "Found mean energy $(round(result.µ, digits=2)) for channel $channel / det $det"
    fig_peakfit = LegendMakie.lplot(report, figsize = (600, 430), titlesize = 17, title = get_plottitle(filekey, det, "Energy Distribution"), juleana_logo = juleana_logo, xlabel = "Energy (ADC)")
    pname = plt_folder * _get_pltfilename(data, filekey, channel, Symbol("peakfit_$(e_type)"))
    save(pname, fig_peakfit )
    @info "Save peak fit plot to $pname"

    # pulser
    e_pulser = filter!(isfinite, dsp_pars.pulser_e_trap)
    qmin, qmax = 0.0, 1.0
    emin = quantile(e_pulser, qmin)
    emax = quantile(e_pulser, qmax)
    cuts_τ = cut_single_peak(e_pulser, emin, emax,; n_bins=nbins, relative_cut=rel_cut_fit)
    result_pulser, report_pulser = fit_single_trunc_gauss(e_pulser, cuts_τ)
    @debug "Found mean energy $(round(result.µ, digits=2)) for channel $channel / det $det"
    fig_pulser = LegendMakie.lplot(report_pulser, figsize = (600, 430), titlesize = 17, title = get_plottitle(filekey, det, "Energy Distribution Pulser"), juleana_logo = false, xlabel = "Energy (ADC)")
    pname_pulser = plt_folder * _get_pltfilename(data, filekey, channel, Symbol("peakfit_$(e_type)_pulser"))
    save(pname_pulser, fig_pulser )
    @info "Save peak fit plot to $pname_pulser"  

    result_ecal = (µ = result.μ, µ_pulser = result_pulser.μ, fit = result, fit_pulser = result_pulser)
    writelprops(data.par[category].rpars.ecal[period], run, PropDict("$channel" => result_ecal))
    @info "Saved pars to disk"
end