function process_noisesweep end
export process_noisesweep
function process_noisesweep(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig; reprocess::Bool = false, filter_type::Symbol = :trap, waveform_type::Symbol = :waveform, n_evts::Union{<:Float64, <:Int} = NaN, scoperun::Bool = false)
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    
    # load or calculate noise sweep 
    pars_pd = if isfile(joinpath(data_path(data.par[category].rpars.noise[period]), "$run.json" )) && !reprocess
        @info "Load noise sweep results for $category-$period-$run-$channel"
        data.par[category].rpars.noise[period,run, channel]
    else 
        PropDict()
    end 
    if !haskey(pars_pd, Symbol(waveform_type))
        # load data and do noise sweep 
        data_raw = read_ldata(data, DataTier(:raw), filekeys, channel)
        wvfs = if isfinite(n_evts)
            n_evts = n_evts < length(data_raw[waveform_type]) ? n_evts : length(data_raw[waveform_type])
            data_raw[waveform_type][1:n_evts]
        else
            data_raw[waveform_type]
        end
        result_rt, report_rt = noise_sweep(filter_type, wvfs, dsp_config)
        result_rt = merge(result_rt, (rt = report_rt.rt, noise = report_rt.noise, ft = report_rt.ft))

        # merge with previous results and save. 
        pars_pd = merge(pars_pd, Dict(Symbol("$(waveform_type)") => result_rt))
        writelprops(data.par[category].rpars.noise[period], run, PropDict("$channel" => pars_pd))
        @info "Save noise sweep results to pars (type $waveform_type, $(filter_type))"
    else
        # read result, and generate report (redo interpolation)
        result_rt = data.par[category].rpars.noise[period, run, channel][waveform_type]
        f_interp = let enc = result_rt.noise, rt = ustrip.(result_rt.rt)
            if length(enc) >= 4
                BSinterpolate(rt, enc, BSplineOrder(4))
            else
                LinearInterpolation(rt, enc)
            end
        end

        report_rt = merge(result_rt, PropDict(:f_interp => f_interp))
    end 

    # plot and save 
    asic_meta = data.metadata.hardware.asic(filekeys[1])

    plot_fun, plot_units = if scoperun
        plot_noise_sweep_osci, [:e, :keV, :µV]
    else
        plot_noise_sweep, [:ADC, :e, :keV, :µV]
    end

    for yunit in plot_units
        plt = plot_fun(report_rt, yunit; 
                DAQ_bits = 14, DAQ_dynamicrange_V  = 2.0, 
                gain =  asic_meta.gain_tot, 
                cap_inj = ustrip.(uconvert(u"F", asic_meta.cap_inj)),
                title = get_plottitle(filekeys[1], _channel2detector(data, channel), "Noise sweep") * @sprintf("\n gain add. = %.1f, Cf = %.0f fF, Cinj = %.0f fF",  asic_meta.gain_tot, ustrip.(uconvert(u"fF", asic_meta.cap_feedback)), ustrip.(uconvert(u"fF", asic_meta.cap_inj)) )) 
        plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], :noise_sweep) * "/"
        pname = plt_folder *  _get_pltfilename(data, filekeys[1], channel, Symbol("noisesweep_$(filter_type)_$(waveform_type)_$(plt.yunit)"))
        save(pname, plt.fig)
        @info "Save plot to $pname"
    end 

    return result_rt, report_rt
end


function plot_noise_sweep(report, yunit::Symbol; DAQ_bits::Int = 14, DAQ_dynamicrange_V::T = 2.0, gain::T = 1.0,  cap_inj::T = 500.0*1e-15, title = "") where T<:Real
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
        title = title,
        xlabel = "Rise time ($x_unit)", ylabel = "Noise$ylabel_extra ($yunit)",
        limits = ((extrema(x)[1] - 0.2, extrema(x)[2] + 0.2), (nothing, nothing))) 
    lines!(ax, x_inter, y_inter, color = :deepskyblue2, linewidth = 3, linestyle = :solid, label =  @sprintf("noise min = %.1f %s (ft = %.2f %s, rt = %.1f %s)", report.min_noise * y_scale, yunit, ustrip(report.ft), unit(report.ft), ustrip(report.rt_opt), unit(report.rt_opt)))
    Makie.scatter!(ax, x, y,  color = :black, label = "Data")
    axislegend()
    fig
    return (fig = fig, ax = ax, y_scale = y_scale, yunit = yunit)
end
export plot_noise_sweep

function plot_noise_sweep_osci(report, yunit::Symbol; gain::T = 1.0,  cap_inj::T = 500.0*1e-15, title = "", DAQ_bits::Int = 1, DAQ_dynamicrange_V::T = 1.0) where T<:Real
    ylabel_extra = "" 
    y_scale = if yunit == :µV
        ylabel_extra = " / gain" ;
        1e6 ./ gain
    elseif yunit == :keV
        V_to_electrons(1.0, cap_inj; gain = gain) *  Ge_Energy_per_eholePair(90) / 1e3
    elseif yunit == :e
        V_to_electrons(1.0, cap_inj; gain = gain)
    else
        @warn("Invalid yunit")
        return
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
        title = title,
        xlabel = "Rise time ($x_unit)", ylabel = "Noise$ylabel_extra ($yunit)",
        limits = ((extrema(x)[1] - 0.2, extrema(x)[2] + 0.2), (nothing, nothing))) 
    lines!(ax, x_inter, y_inter, color = :deepskyblue2, linewidth = 3, linestyle = :solid, label =  @sprintf("noise min = %.1f %s (ft = %.2f %s, rt = %.1f %s)", report.min_noise * y_scale, yunit, ustrip(report.ft), unit(report.ft), ustrip(report.rt_opt), unit(report.rt_opt)))
    Makie.scatter!(ax, x, y,  color = :black, label = "Data")
    axislegend()
    fig
    return (fig = fig, ax = ax, y_scale = y_scale, yunit = yunit)
end
export plot_noise_sweep_osci