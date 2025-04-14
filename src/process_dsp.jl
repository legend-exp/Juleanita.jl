
"""
    process_dsp(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, τ_pz::Quantity{<:Real}, pars_filter::PropDict; reprocess::Bool = false )  
    process_dsp(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId ; kwargs... ) 

Goal: 
- run the DSP processing for all raw files in the given period, run, category and channel.
- based on `simple_dsp` function
- save the results in the jldsp tier
- if reprocess is false, it will skip the files that are already processed

INPUTS:
    - `data::LegendData` LegendData object. You need `"LEGEND_DATA_CONFIG"` to construct this, e.g. `l200 = LegendData(:l200)`
    - `period::DataPeriod` data period, e.g. `DataPeriod(1)`
    - `run::DataRun` data run, e.g. `DataRun(1)`
    - `category::Symbol` data category, e.g. `DataCategory(:cal)`
    - `channel::ChannelId` channel id, e.g. `ChannelId(1)` (depending on your data!)
    - `dsp_config::DSPConfig` DSP configuration object. If not specified will take default from metadata
    - `τ_pz::Quantity{<:Real}` decay time used for pole-zero correction. If not specified will take from rpars.pz
    - `pars_filter::PropDict` optimized filter parameters used in DSP. If not specified will take from rpars.fltopt

KWARGS:
    - `reprocess::Bool` reprocess the files or not
    
OUTPUTS:
    - save the DSP results in the jldsp tier
    - print the progress
    - print the completion message
"""
function process_dsp end
export process_dsp
function process_dsp(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, τ_pz::Quantity{<:Real}, pars_filter::PropDict; reprocess::Bool = false)  

    @debug "Start DSP processing for period $period, run $run, channel $channel"
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    dsp_folder =  data.tier[DataTier(:jldsp), category , period, run] * "/"
    if reprocess == false
        dsp_files = dsp_folder .* string.(filekeys) .* "-tier_jldsp.lh5"
        if sum(isfile.(dsp_files)) == length(filekeys)
            @info "all raw files are already processed. You're already done!"
            println("to read dsp files: use read_ldata(data, :jldsp, category , period, run, channel)")
            return 
        elseif sum(isfile.(dsp_files)) > 0 
            @info "$(sum(isfile.(dsp_files))) raw files are already processed - skip these files"
        end
        filekeys = filekeys[.!isfile.(dsp_files)]
    end

     # save dsp results
    if !ispath(dsp_folder)
        mkpath(dsp_folder)
    end

   # Threads.@threads 
    for f in eachindex(filekeys) 
        filekey = filekeys[f]
        data_raw = read_ldata(data, DataTier(:raw), filekey, channel)
        @debug "processing waveform data"
        dsp_par = simple_dsp(Table(data_raw), dsp_config; τ_pz = τ_pz, pars_filter = pars_filter)
       
        if hasproperty(data_raw, :pulser)
            @debug "processing pulser data"
            dsp_par_pulser = simple_dsp_pulser(Table(data_raw), dsp_config; τ_pz = 0.0u"µs", pars_filter = pars_filter)
        end

        dsp_file = dsp_folder * string(filekey) * "-tier_jldsp.lh5"
        fdsp = lh5open(dsp_file, "w")
        for par in columnnames(dsp_par)
            fdsp["$channel/jldsp/$par"]  = getproperty(dsp_par, par)
        end
        if hasproperty(data_raw, :pulser)
            for par in columnnames(dsp_par_pulser)
                fdsp["$channel/jldsp/$par"]  = getproperty(dsp_par_pulser, par)
            end
        end
        fdsp["$channel/jldsp/eventnumber"] = data_raw.eventnumber
        fdsp["$channel/jldsp/timestamp"] = data_raw.timestamp

        close(fdsp)
        println("dsp processing done for $filekey")
    end
    println("DONE! To read the generated dsp files: read_ldata(data, :jldsp, category , period, run, channel)")
end

function process_dsp(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId ; kwargs... )  
    @info "use default DSP config and filter parameter "
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    dsp_config = DSPConfig(dataprod_config(data).dsp(filekeys[1]).default)
    τ_pz = mvalue(get_values(data.par[category].rpars.pz[period, run, channel]).τ)
    pars_filter = data.par[category].rpars.fltopt[period,run,channel]
    process_dsp(data, period, run, category, channel, dsp_config, τ_pz, pars_filter; kwargs...)
end 