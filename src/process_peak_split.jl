"""
    process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, ecal_config::PropDict, dsp_config::DSPConfig, qc_config::PropDict; reprocess::Bool = false)
    process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
Create peak files containting only waveforms from peaks in the calibration spectrum.
1. Read raw data
2. Find peaks in the calibration spectrum using rough energy estimate `data_fk.daqenergy`
3. Do simple DSP for peaks only 
4. apply quality cuts based on simple DSP
5. Save waveforms after QC to peak files
## inputs:
- data: LegendData object
- period: DataPeriod object
- run: DataRun object
- category: Symbol or DataCategory object, e.g. :cal for calibration 
- channel: ChannelId for germanium detector 
- ecal_config: PropDict, energy calibration configuration
- dsp_config: DSPConfig, DSP configuration
- qc_config: PropDict, quality cut configuration

## keyword arguments:
reprocess: Bool, default is false
## Output:
- h_uncals histograms of peaksearch
- peakpos 
Also, peak files containting only waveforms from peaks in the calibration spectrum are saved to jlpeaks folder.
"""
function process_peak_split end
export process_peak_split

function process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, ecal_config::PropDict, dsp_config::DSPConfig, qc_config::PropDict; reprocess::Bool = false, plotHist::Bool = true, peakfinder_gamma::Union{Int, Float64} = NaN)
                        
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    peak_folder =  data.tier[DataTier(:jlpeaks), category , period, run] * "/"
    if !ispath(peak_folder)
        mkpath(peak_folder)
        @info "create path: $peak_folder"
    end
    peak_files = peak_folder .* string.(filekeys) .* "-tier_jlpeaks.lh5"

    if reprocess == false
        if sum(isfile.(peak_files)) == length(filekeys)
            @info "all peak files are already processed. You're already done!"
            println("to read peak files: use read_ldata(data, :jlpeaks, category , period, run, channel)")
            return 
        else
            @info "$(sum(isfile.(peak_files))) peak files are already processed - skip"
        end
        filekeys = filekeys[.!isfile.(peak_files)]
    end

    if Symbol(ecal_config.source) == :co60
        gamma_lines =  ecal_config.co60_lines
        gamma_names =  ecal_config.co60_names
        left_window_sizes = ecal_config.co60_left_window_sizes 
        right_window_sizes = ecal_config.co60_right_window_sizes 
    elseif Symbol(ecal_config.source) == :th228
        gamma_lines =  ecal_config.th228_lines[end]
        gamma_names =  [ecal_config.th228_names[end]]
        left_window_sizes = 1.5 * ecal_config.th228_left_window_sizes[end]
        right_window_sizes = 1.5 * ecal_config.th228_right_window_sizes[end] 
    end

    if peakfinder_gamma !== NaN
        gamma_lines =  gamma_lines[peakfinder_gamma]
        gamma_names =  gamma_names[peakfinder_gamma]
        left_window_sizes = left_window_sizes[peakfinder_gamma]
        right_window_sizes = right_window_sizes[peakfinder_gamma]
        @info "use only peakfinder_gamma $(peakfinder_gamma) for peak search"
    end

    # result_peaksearch = Dict()
    function _search_peaks(X::Vector{<:Real}; peaks::Union{<:Quantity, Vector{<:Quantity}} = gamma_lines)
        # binning and peak search windows and histogram settings 
        bin_min = quantile(X, ecal_config.left_bin_quantile)
        bin_max = quantile(X, ecal_config.right_bin_quantile)
        peak_min = quantile(X, ecal_config.left_peak_quantile)
        peak_max = quantile(X, ecal_config.right_peak_quantile)
        bin_width = LegendSpecFits.get_friedman_diaconis_bin_width(filter(in(bin_min..bin_max), X))
        if (peak_max-peak_min)/bin_width < ecal_config.nbins_min
            bin_width = (peak_max-peak_min)/ ecal_config.nbins_min
        end

        # peak search
        peakpos = []
        function _peakfinder(p_min::T, p_max::T) where T <: Real
            try 
                h_peaksearch = StatsBase.fit(Histogram, e_uncal, p_min:bin_width:p_max) # histogram for peak search
                _, peakpos = RadiationSpectra.peakfinder(h_peaksearch, σ= ecal_config.peakfinder_σ, backgroundRemove=true, threshold = ecal_config.peakfinder_threshold)
                @info "Found $(length(peakpos)) peak(s) at $(peakpos)"
                return peakpos
            catch e
                @info "peakfinder failed: $e"
                return NaN
            end
        end
        
        peakpos = _peakfinder(peak_min, peak_max)
        if peakpos == NaN
            @debug "try larger peak window"
            peakpos = _peakfinder(peak_min, 1.5*peak_max)  
        elseif length(peakpos) > length(peaks) 
            @debug"too many peaks found - try smaller peak window"
            for i in 0.1:0.1:0.9
                peakpos = _peakfinder((1+i)*peak_min, (1-i)*peak_max)
                if length(peakpos) == length(peaks) 
                    break 
                end 
            end 
        elseif length(peakpos) < length(peaks)
          @debug "too few peaks found - try larger peak window"
            for i in 0.1:0.1:0.9
                peakpos = _peakfinder((1-i)*peak_min, (1+i)*peak_max)
                if length(peakpos) == length(peaks) 
                    break 
                end 
            end 
        end 

        if (length(peakpos) !== length(peaks)) || peakpos == NaN
            error("Number of peaks found $(length(peakpos)) ($peakpos); expected gamma lines $(length(peaks)) \n you could try to modify peakfinder_threshold and/or peakfinder_σ")
        end 
        cal_simple = mean(peaks./sort(peakpos))
        e_cal = X .* cal_simple 
        result = (e_simplecal = e_cal, peakpos = peakpos, hist_bins = 0:bin_width:maximum(X), cal_simple = cal_simple)
        return result 
    end

    # do peaksearch on all filekeys 
    e_uncal =  filter(x-> x >= qc_config.e_trap.min, read_ldata(:daqenergy, data, DataTier(:raw), filekeys, channel));
    if isempty(e_uncal)
         @error "No energy values >= $(qc_config.e_trap.min) found"
    end
    result_ps =  _search_peaks(e_uncal; peaks = gamma_lines);  # do peak search
    @debug "peak search done"
    
    # apply simple calibration to all filekeys  
   if plotHist == true 
        plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], :peak_split) * "/"
        if !ispath(plt_folder)
            mkpath(plt_folder)
            @debug "create path: $plt_folder"
        end

        if period >= DataPeriod(3)
            xunit = "ADC"
            xformatter = X -> map(x -> "$(Int(x))", X)
        else
            xunit = "V"
            xformatter = X -> map(x -> "$(round(x, digits = 2))", X)
        end

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel = "Energy ($xunit)", ylabel = "Counts", xtickformat = xformatter, title = get_plottitle(filekeys[1], _channel2detector(data, channel), "peak split"), limits = ((nothing, nothing), (0, nothing)))
        Makie.stephist!(ax,result_ps.e_simplecal ./ result_ps.cal_simple, bins = result_ps.hist_bins  )
        Makie.hist!(ax, result_ps.e_simplecal ./ result_ps.cal_simple, bins = result_ps.hist_bins) 
        vlines!(ax, result_ps.peakpos, color = :red2, label =  "peakfinder result", alpha = 0.7, linestyle = :dash, linewidth = 2.0)
        axislegend()
        fig  
        pname = plt_folder * "peak_split.png"
        save(pname, fig)
    end  

    for f in eachindex(filekeys) 
        filekey = filekeys[f]
        peak_file = peak_files[f]

        # read raw data (waveform tier) from file 
        data_fk = read_ldata(data, DataTier(:raw), filekey, channel)
        e_uncal = filter(x -> x >= qc_config.e_trap.min , data_fk.daqenergy)
        if isempty(e_uncal)
            @warn "No energy values >= $(qc_config.e_trap.min) found for $filekey - skip"
            continue
        end

        # # do peak search
        # result_ps =  _search_peaks(e_uncal; peaks = gamma_lines);
        e_simplecal = e_uncal .* result_ps.cal_simple

        # save results to peakfile 
        wvfs = data_fk.waveform
        eventnumber = data_fk.eventnumber
        timestamp = data_fk.timestamp
        
        # save
       
        fid = lh5open(peak_file, "w") 
        for i in eachindex(gamma_lines)
            peakIdx = findall((gamma_lines[i] - left_window_sizes[i]) .<= e_simplecal .< (gamma_lines[i] + right_window_sizes[i]))
            # do simple dsp. only for peaks. 
            dsp_par = simple_dsp_qc(Table(waveform = wvfs[peakIdx]), dsp_config)
            qc_cuts = apply_qc(dsp_par, qc_config)
            qc_flag = qc_cuts.wvf_keep.all
            qc_idx = findall(qc_flag)
            qc_surv = qc_cuts.qc_surv.all
            fid["$channel/jlpeaks/$(gamma_names[i])/waveform"] = wvfs[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/daqenergy"] = e_uncal[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/e_simplecal"] = e_simplecal[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/eventnumber"] = eventnumber[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/timestamp"] = timestamp[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/gamma_line"] = fill(gamma_lines[i], length(qc_idx))
            fid["$channel/jlpeaks/$(gamma_names[i])/qc_flag"] = qc_flag[qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/qc_surv"] = fill(qc_surv, length(qc_idx))
            for par in columnnames(dsp_par)
                fid["$channel/jlpeaks/$(gamma_names[i])/$par"]  = getproperty(dsp_par, par)
            end
        end
        println("peak file processing done for $filekey")
        close(fid)
    end

end

function process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    qc_config = dataprod_config(data).qc(filekeys[1]).default
    ecal_config = dataprod_config(data).energy(filekeys[1]).default
    dsp_config = DSPConfig(dataprod_config(data).dsp(filekeys[1]).default)
    @info "use default configs: qc, ecal, dsp"
    process_peak_split(data, period, run, category, channel, ecal_config, dsp_config, qc_config; kwargs...)
end

