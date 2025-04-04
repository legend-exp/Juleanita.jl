

"""
    process_ctc(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, ecal_config::PropDict, ctc_config::PropDict;
                energy_types::Vector{Symbol} = Symbol.(ctc_config.energy_types), reprocess::Bool = false)
    process_ctc(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
Calculate charge-trapping correction: by looking at correlation between drift time and energy. Save correction function to rpars. 
Inputs:
- `data::LegendData`: LegendData object
- `period::DataPeriod`: data period
- `run::DataRun`: data run
- `category::Union{Symbol, DataCategory}`: data category, e.g. :cal
- `channel::ChannelId`: channel id
- `ecal_config::PropDict`: energy calibration configuration. If not specified use default from metadata
- `ctc_config::PropDict`: charge-trapping correction configuration. If not specified use default from metadata
Optional:
- `energy_types::Vector{Symbol}`: energy types to process
- `reprocess::Bool`: reprocess the files or not
- `juleana_logo::Bool`: add juleana logo to plots
"""
function process_ctc(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, ecal_config::PropDict, ctc_config::PropDict; 
   energy_types::Vector{Symbol} = Symbol.(ctc_config.energy_types), juleana_logo::Bool = false, reprocess::Bool = false)

    @debug "Create pars db"
    mkpath(joinpath(data_path(data.par[category].rpars.ctc), string(period)))
    pars_db = PropDict(data.par[category].rpars.ctc[period, run])
    pars_db = ifelse(reprocess, PropDict(), pars_db)
    if reprocess @info "Start reprocessing" end

    result_dict    = Dict{Symbol, NamedTuple}()
    processed_dict = Dict{Symbol, Bool}()
    plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekey, :ctc) * "/"

    # check if energy types are already processed
    if !reprocess && haskey(pars_db, det)
        @debug "Channel $(det) already processed, check missing energy types"
        for e_type in energy_types
            if haskey(pars_db[det], e_type)
                @debug "Filter $e_type already processed, skip"
                log_info = log_nt((channel, det, ProcessStatus(1), e_type, pars_db[det][e_type].fct*1e6, pars_db[det][e_type].fwhm_before, pars_db[det][e_type].fwhm_after, "Already processed --> skipped."))
                processed_dict[e_type] = false
                log_info_dict[e_type] = log_info
            end
        end
    end

    # get simple calibration config
    source = Symbol(ecal_config.source)
    @debug "Load simple calibration configuration for $(source) source"
    if source == :th228
        calib_type = :th228
    elseif source == :co60
        calib_type = :gamma
    else
        error("Unknown source $(source)")
    end
    gamma_lines =  ecal_config[Symbol("$(source)_lines")]
    left_window_sizes = ecal_config[Symbol("$(source)_left_window_sizes")]
    right_window_sizes = ecal_config[Symbol("$(source)_right_window_sizes")]

    # start processing
    for e_type in energy_types
        if haskey(processed_dict, e_type)
            continue
        end
        
        # load data: uncalibrated energies and qdrift
        e_uncal = nothing
        qdrift = nothing
        try
            @debug "Load uncalibrated energies after qc"
            dsp_par = read_ldata(data, :jldsp, category, period, run, channel)
            e_uncal = getproperty(dsp_par, e_type)[dsp_par.qc]
            qdrift = getproperty(dsp_par, :qdrift)[dsp_par.qc]
        catch e
            @error "Error in loading data for channel $channel: $e"
            throw(ErrorException("Error data loader"))
        end
     
        try
            @debug "Correct $e_type"
            # do simple calibration
            result_simple, report_simple = nothing, nothing
            try
                @debug "Get $e_type simple calibration"
                result_simple, report_simple = simple_calibration(e_uncal, gamma_lines , left_window_sizes, right_window_sizes,; 
                calib_type = calib_type, binning_peak_window=ecal_config.binning_peak_window, quantile_perc=NaN, 
                peak_quantile= ecal_config.left_peak_quantile..ecal_config.right_peak_quantile, 
                bin_quantile = ecal_config.left_bin_quantile..ecal_config.right_bin_quantile, 
                peakfinder_threshold = ecal_config.peakfinder_threshold, 
                peakfinder_σ = ecal_config.peakfinder_σ);
            catch e
                @error "Error in $e_type simple calibration for channel $channel: $(e)"
                throw(ErrorException("Error in $e_type simple calibration"))
            end
            # get simple calibration constant
            m_cal_simple = result_simple.c
            # save plots for simple calibration for control
            report_simple_alt = (h_calsimple = report_simple.h_calsimple, 
                h_uncal = report_simple.h_uncal,
                c = report_simple.c,
                fep_guess = report_simple.peak_guess,
                peakhists = report_simple.peakhists,
                peakstats = report_simple.peakstats)
            fig_simple = LegendMakie.lplot(report_simple_alt, title = get_plottitle(filekey, det, "Simple Calibration"; additiional_type=string(e_type)), cal = true, juleana_logo = juleana_logo)
            Makie.current_axis().titlesize = 16;
            [delete!(leg) for leg in fig_simple.content if leg isa Legend]
            vl = vlines!([report_simple.peak_guess * ustrip(report_simple.c)], color = :red2, label = "Peak Guess", alpha = 0.5, linewidth = 3)
            axislegend(Makie.current_axis(), [vl], ["Peak Guess"], position = :lt)
            Makie.xlims!(Makie.current_axis(), 300, 2000)
            pname_simple = plt_folder * _get_pltfilename(data, filekey, channel, Symbol("simple_calibration_ctc_$(e_type)"))
            save(pname_simple, fig_simple)

            # do charge-trapping correction 
            result_ctc, report_ctc = nothing, nothing
            try
                @debug "Get $e_type Charge Trapping Alpha"
                result_ctc, report_ctc = ctc_energy(e_uncal .* m_cal_simple, qdrift, ctc_config.peak, (ctc_config.left_window_size, ctc_config.right_window_size), m_cal_simple; e_expression="$e_type", pol_order=ctc_config.ctc_order)
            catch e
                @error "Error in $e_type alpha generation $channel: $e"
                throw(ErrorException("Error in $e_type alpha generation"))
            end
            @debug "Found Best $e_type FWHM: $(round(u"keV", result_ctc.fwhm_after, digits=2))"
            @debug "Found $e_type FCTs: $(round.(result_ctc.fct .* 1e6, digits=2))E-6"
            
            fig_ctc = LegendMakie.lplot(report_ctc, figsize = (600,600), title = get_plottitle(filekey, det, "CTC"; additiional_type="$e_type $(ctc_config.peak)"), juleana_logo = juleana_logo, watermark = false)
            Makie.current_axis().titlesize = 16;
            Makie.ylims!(Makie.current_axis(), 0, quantile(qdrift./maximum(qdrift), 0.995))
            fig_ctc
            pname_ctc = plt_folder * _get_pltfilename(data, filekey, channel, Symbol("ctc_$(e_type)"))
            save(pname_ctc, fig_ctc)

            # add results to dict
            result_dict[e_type]   = result_ctc
            processed_dict[e_type] = true
        catch e
            @error "Error in $e_type CT correction: $e"
            # add results to dict
            processed_dict[e_type] = false
        end
    end

    @info "Finished CT correction"

    # save results to pars 
    result_ctc = PropDict(Dict("$channel" => result_dict))
    writelprops(data.par[category].rpars.ctc[period], run, result_ctc)
    writevalidity(data.par[category].rpars.ctc, filekey, (period, run))
    @info "Saved pars to disk"
    return (result = result_ctc, status = processed_dict)
end 

function process_ctc(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
    @info "use default ecal and ctc configs"
    filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
    ecal_config = dataprod_config(data).energy(filekey).default
    ctc_config = dataprod_config(data).energy(filekey).ctc.default
    process_ctc(data, period, run, category, channel, ecal_config, ctc_config; kwargs...)
end
 