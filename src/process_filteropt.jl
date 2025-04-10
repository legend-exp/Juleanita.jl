"""
    process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, τ_pz::Quantity{T}, peak::Symbol; rt_opt_mode::Symbol = :blnoise, reprocess::Bool = false, filter_types::Vector{Symbol} = [:trap, cusp])
    process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)            

Filter optimization for filter_types
- load waveforms from peakfile, shift baseline and pole-zero 
- optimize rise-time for minimum baseline noise after filtering
- optimize flat-top time for FWHM of peak
- save results to disk
- sanity plots for rise-time and flat-top time optimization
Inputs: 
- `data::LegendData`: LegendData object
- `period::DataPeriod`: data period
- `run::DataRun`: data run
- `category::Union{Symbol, DataCategory}`: data category, e.g. :cal
- `channel::ChannelId`: channel id
Optional:
- `dsp_config::DSPConfig`: DSP configuration object. If not specified will take default from metadata
- `τ_pz::Quantity{T}`: pole-zero decay time.  If not specified will take default from metadata
- `peak::Symbol`: peak to optimize for (needs existing peakfile!). If not specified will take default from metadata
- `rt_opt_mode::Symbol`: mode for rise-time optimization (:blnoise or :pickoff) --> two different strategies for optimization
- `reprocess::Bool`: reprocess the files or not
- `filter_types::Vector{Symbol}`: filter types to optimize for
"""
function process_filteropt end
export process_filteropt
function process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, τ_pz::Quantity{T}, peak::Symbol; 
                rt_opt_mode::Symbol = :bl_noise, reprocess::Bool = false, filter_types::Vector{Symbol} = [:trap, :cusp, :zac], fwhm_rel_cut_fit::T = 0.1) where T<:Real 
    det = _channel2detector(data, channel)
    @info "Optimize filter for period $period, run $run, channel $channel /det $det - $filter_types"

    # check if decaytime pars already exist
    fltopt_file = joinpath(mkpath(data_path(data.par[category].rpars.fltopt[period])), "$(string(run)).json")
    if isfile(fltopt_file) && !reprocess
        @info "Filter optimization file already exist for $category period $period - run $run - channel $channel - you're done!"
        return
    end

    # prepare results dict
    mkpath(joinpath(data_path(data.par[category].rpars.fltopt), string(period)))
    result_filteropt_dict = Dict{Symbol, NamedTuple}()
    @debug "Created path for filter optimization results"

    # load waveforms from peakfile
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    if peak == :all 
        data_peak = read_ldata(data, DataTier(:raw), filekeys, channel)
        data_peak = if Symbol(category) == :cal 
            merge(data_peak, (gamma_line = [1170*u"keV"],))
        elseif Symbol(category) == :bch
            merge(data_peak, (gamma_line = [1000*u"keV"],))
        end
    else
        data_peak  = read_ldata((peak), data, :jlpeaks, category, period, run, channel)
    end 
    wvfs = data_peak.waveform
    @debug "Loaded waveforms for peak $peak"

    function process_filteropt_fltr(filter_type::Symbol)
        @info "Optimize filter $filter_type"
        _, def_ft = get_fltpars(PropDict(), filter_type, dsp_config)
        filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
        plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekey, :filteropt) * "/"

        # STEP 1: rise-time optimization --> min. baseline noise after filtering 
        if rt_opt_mode == :bl_noise 
            result_rt, report_rt = filteropt_rt_optimization_blnoise(filter_type, wvfs, dsp_config, τ_pz; ft = def_ft)
           
            rt_inter = range(ustrip.(report_rt.rt[1]), stop = ustrip(maximum(report_rt.rt[findall(isfinite.(report_rt.noise))])), step = 0.05); 
            p = Figure()
            ax = Axis(p[1, 1], 
                xlabel = "Rise time ($(unit(report_rt.rt_opt)))", ylabel = "Noise (a.u.)",
                limits = ((ustrip.(extrema(report_rt.rt))[1] - 0.2, ustrip.(extrema(report_rt.rt))[2] + 0.2), (nothing, nothing)),
                title = "Noise sweep ($filter_type), $period-$run-$channel, $peak peak \n" * @sprintf("fixed ft = %.2f %s, optimal rt = %.1f %s", ustrip(def_ft), unit(def_ft), ustrip(report_rt.rt_opt), unit(report_rt.rt_opt)), )
            lines!(ax, rt_inter, report_rt.f_interp.(rt_inter), color = :deepskyblue2, linewidth = 3, linestyle = :solid, label = "Interpolation")
            Makie.scatter!(ax, ustrip.(collect(report_rt.rt)), report_rt.noise,  color = :black, label = "Data")
            axislegend()
            pname = plt_folder * split(LegendDataManagement.LDMUtils.get_pltfilename(data, filekeys[1], channel, Symbol("noise_sweep_$(filter_type)_blnoise")),"/")[end]
            d = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], Symbol("noise_sweep_$(filter_type)_blnoise"))
            ifelse(isempty(readdir(d)), rm(d), nothing )
        elseif rt_opt_mode == :pickoff
            # for now: do 2 versions or rt...the same using LEGEND DSP function 
            # pro: + compatibility with Juleana; proper fitting of enc instead of rms
            # con: - doesnt work well for smaller number of waveforms. result seems to be more noisy, slower?
            enc_grid = getfield(Main, Symbol("dsp_$(filter_type)_rt_optimization"))(wvfs, dsp_config,  τ_pz; ft=def_ft)
            enc_min, enc_max = _quantile_truncfit(enc_grid; qmin = 0.02, qmax = 0.98)
            e_grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
            result_rt, report_rt = fit_enc_sigmas(enc_grid, e_grid_rt, enc_min, enc_max, round(Int,size(enc_grid)[2]/5), 0.1)
            @info "Found optimal rise-time: $(result_rt.rt) at fixed ft = $def_ft" 
            p = LegendMakie.lplot(report_rt, title = get_plottitle(filekeys[1], det, "Noise Sweep"; additiional_type=string(filter_type)))
            pname = plt_folder * split(LegendDataManagement.LDMUtils.get_pltfilename(data, filekeys[1], channel, Symbol("noise_sweep_$(filter_type)_pickoff")),"/")[end]
            d = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], Symbol("noise_sweep_$(filter_type)_pickoff"))
            ifelse(isempty(readdir(d)), rm(d), nothing )
        end 
        save(pname, p)
        @info "Save sanity plot to $pname"

        # 2. flat top time optimixation 
        ft_qmin, ft_qmax = dsp_config.kwargs_pars.ft_qmin, dsp_config.kwargs_pars.ft_qmax
        e_grid_ft   = getproperty(dsp_config, Symbol("e_grid_ft_$(filter_type)"))
        e_grid = getfield(LegendDSP, Symbol("dsp_$(filter_type)_ft_optimization"))(wvfs, dsp_config, τ_pz, mvalue(result_rt.rt_opt))
        e_min, e_max = _quantile_truncfit(e_grid; qmin = ft_qmin, qmax = ft_qmax)
        result_ft, report_ft = fit_fwhm_ft(e_grid, e_grid_ft, result_rt.rt_opt,  e_min, e_max, fwhm_rel_cut_fit; peak = data_peak.gamma_line[1])
        @info "Found optimal flattop-time: $(result_ft.ft) with FWHM $(round(u"keV", result_ft.min_fwhm, digits=2))"
        p = Figure()
        LegendMakie.lplot!(report_ft, title = get_plottitle(filekey, det, "$peak FT Scan"; additiional_type=string(filter_type)), juleana_logo = false)
        if Symbol(category) == :bch
            axs =  p.content[findall(map(x -> x isa Axis, p.content))]
            axs[1].ylabel = "FWHM (a.u.)"
        end
        pname_ft = plt_folder * split(LegendDataManagement.LDMUtils.get_pltfilename(data, filekeys[1], channel, Symbol("fwhm_ft_scan_$(filter_type)")),"/")[end]
        d = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], Symbol("fwhm_ft_scan_$(filter_type)"))
        ifelse(isempty(readdir(d)), rm(d), nothing )
        save(pname_ft, p)
        @info "Save sanity plot to $pname_ft"

        # return result for this filter type
        merge(result_rt, result_ft)
    end 

    # filter optimization: rise-time, and flat-top times
    for filter_type in filter_types
        result_filteropt_dict[filter_type] =  process_filteropt_fltr(filter_type)
    end
    result = PropDict(Dict("$channel" => result_filteropt_dict))
    writelprops(data.par[category].rpars.fltopt[period], run, result)
    @info "Saved pars to disk"
end

function process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
    
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    dsp_config = DSPConfig(dataprod_config(data).dsp(filekeys[1]).default)
    pz_config = dataprod_config(data).dsp(filekeys[1]).pz.default

    peak =  Symbol(pz_config.peak)
    τ_pz = mvalue(get_values(data.par[category].rpars.pz[period, run, channel]).τ)
    @debug "Loaded decay time for pole-zero correction: $τ_pz"

    process_filteropt(data, period, run, category, channel, dsp_config, τ_pz, peak; kwargs...) 
end