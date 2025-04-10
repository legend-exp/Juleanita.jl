function process_noisesweep end
export process_noisesweep
function process_noisesweep(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig; reprocess::Bool = false, filter_type::Symbol = :trap, waveform_type::Symbol = :waveform, diff_output::Bool = true, n_evts::Int = 1000)
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    
    # load or calculate noise sweep 
    pars_pd = if isfile(joinpath(data_path(data.par[category].rpars.noise[period]), "$run.json" )) && !reprocess
        @info "Noise sweep already done for $category-$period-$run-$channel $waveform_type - load results!"
        data.par[category].rpars.noise[period,run, channel]
    else 
        PropDict()
    end 
    if !haskey(pars_pd, Symbol(waveform_type))
        # load data and do noise sweep 
        data_raw = read_ldata(data, DataTier(:raw), filekeys, channel)
        wvfs = data_raw[waveform_type][1:n_evts]
        result_rt, _ = noise_sweep(filter_type, wvfs, dsp_config)
    
        # merge with previous results and save. 
        pars_pd = merge(pars_pd, Dict(Symbol("$(waveform_type)") => result_rt))
        writelprops(data.par[category].rpars.noise[period], run, PropDict("$channel" => pars_pd))
        @info "Save noise sweep results to pars"
    end 

    # read result, and redo interpolation
    result_rt = data.par[category].rpars.noise[period, run, channel][waveform_type]
    f_interp = let enc = result_rt.noise, rt = ustrip.(result_rt.rt)
        if length(enc) >= 4
           BSinterpolate(rt, enc, BSplineOrder(4))
        else
            LinearInterpolation(rt, enc)
        end
    end
    result_rt[:f_interp] = f_interp

    # plot and save 
    asic_meta = data.metadata.hardware.asic(filekeys[1])
    add_gain = if diff_output
        2 * asic_meta.add_gain
    else
        asic_meta.add_gain
    end
    for yunit in [:ADC, :e, :keV, :µV]
        plt = plot_noise_sweep(result_rt, yunit; 
                DAQ_bits = 14, DAQ_dynamicrange_V  = 2.0, 
                add_gain = add_gain, 
                cap_feedback =  ustrip.(uconvert(u"F", asic_meta.cap_feedback)), 
                cap_inj = ustrip.(uconvert(u"F", asic_meta.cap_inj))) 
        plt.ax.title = get_plottitle(filekeys[1], _channel2detector(data, channel), "Noise sweep") * @sprintf("\n gain add. = %.1f, Cf = %.0f fF, Cinj = %.0f fF", add_gain, ustrip.(uconvert(u"fF", asic_meta.cap_feedback)), ustrip.(uconvert(u"fF", asic_meta.cap_inj)) )
        plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], :noise_sweep) * "/"
        pname = plt_folder *  _get_pltfilename(data, filekeys[1], channel, Symbol("noisesweep_$(filter_type)_$(waveform_type)_$(plt.yunit)"))
        save(pname, plt.fig)
        @info "Save plot to $pname"
    end 

    return result_rt
end


function plot_noise_sweep(report, yunit::Symbol; DAQ_bits::Int = 14, DAQ_dynamicrange_V::T = 2.0, add_gain::T = 1.0, cap_feedback::T = 500.0*1e-15, cap_inj::T = 500.0*1e-15) where T<:Real
    gain = add_gain * cap_inj / cap_feedback

    ylabel_extra = "" 
    y_scale = if yunit == :ADC 
        1.0
    elseif yunit == :keV
        _ADC_to_keV(1.0, cap_inj; bits = DAQ_bits, dynamicrange_V = DAQ_dynamicrange_V, gain = gain) 
    elseif yunit == :e
        _ADC_to_electrons(1.0, cap_inj; bits = DAQ_bits, dynamicrange_V = DAQ_dynamicrange_V, gain = gain)
    elseif yunit == :µV
        ylabel_extra = " / gain" ;
        1e6 .* _ADC_to_V(1.0, DAQ_dynamicrange_V, DAQ_bits) ./ gain
    else
        error("Invalid yunit")
    end 
    x = report.rt
    x_unit = unit(x[1])
    x = ustrip.(x)
    y = report.noise .* y_scale
    x_inter = range(x[1], stop = maximum(x[findall(isfinite.(y))]), step = 0.05); 
    y_inter = report.f_interp.(x_inter) .* y_scale
 
    # plot  
    fig = Figure()
    ax = Axis(fig[1, 1], 
        xlabel = "Rise time ($x_unit)", ylabel = "Noise$ylabel_extra ($yunit)",
        limits = ((extrema(x)[1] - 0.2, extrema(x)[2] + 0.2), (nothing, nothing))) 
    lines!(ax, x_inter, y_inter, color = :deepskyblue2, linewidth = 3, linestyle = :solid, label =  @sprintf("noise min = %.1f %s (ft = %.2f %s, rt = %.1f %s)", report.min_noise * y_scale, yunit, ustrip(report.ft), unit(report.ft), ustrip(report.rt_opt), unit(report.rt_opt)))
    Makie.scatter!(ax, x, y,  color = :black, label = "Data")
    axislegend()
    return (fig = fig, ax = ax, y_scale = y_scale, yunit = yunit)
end
export plot_noise_sweep