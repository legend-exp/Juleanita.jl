"""
    processing_hit(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; reprocess::Bool = false, e_types::Vector{<:Symbol} = [:e_trap, :e_cusp, :e_zac])
- run the hit processing for all dsp files in the given period, run, category and channel.
- apply energy calibration to the dsp energy estimator and save the results in the jlhit tier 
- no PSD at the moment, will be added in future 
- if reprocess is false, it will skip the files that are already processed
INPUTS:
- `data::LegendData` LegendData object
- `period::DataPeriod` data period
- `run::DataRun` data run
- `category::Union{Symbol, DataCategory}` data category, e.g. :cal
- `channel::ChannelId` channel id
- `reprocess::Bool` reprocess the files or not
- `e_types::Vector{<:Symbol}` energy types to process. default is [:e_trap, :e_cusp, :e_zac]
OUTPUTS:
- save the hit results in the jlhit tier
"""
function process_hit(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; reprocess::Bool = false, e_types::Vector{<:Symbol} = [:e_trap, :e_cusp, :e_zac])
    hit_folder =  data.tier[DataTier(:jlhit), category , period, run] * "/"
    filekeys = search_disk(FileKey, data.tier[:jldsp, category , period, run])
    hit_files = hit_folder .* string.(filekeys) .* "-tier_jlhit.lh5"
    
    if !reprocess && all(isfile.(hit_files))
        @info "all hit files are already processed. You're done!"
        println("To read the hit files: read_ldata(data, :jlhit, category, period, run, channel)")
        return
    end 

    # load energy calibration
    ecal_pd = data.par[category].rpars.ecal[period, run, channel]

    for f in eachindex(filekeys)
        filekey = filekeys[f]
        hit_file = hit_files[f]

        if !reprocess && isfile(hit_file)
            @info "hit file for filekey $filekey is already processed - skip "
            continue
        end 

        # load dsp parameters for given filekey 
        dsp_pars = Table(read_ldata(data, :jldsp, filekey, channel);)

        function apply_cal(e_type::Symbol)
            calib_func = ljl_propfunc(ecal_pd[e_type].cal.func)
            calib_func.(dsp_pars)
        end

        # save hit results
        if !ispath(hit_folder)
            mkpath(hit_folder)
        end
      
        fhit = lh5open(hit_file, "w")
        for e_type in e_types
            fhit["$channel/jlhit/$e_type"]  = apply_cal(e_type)
        end
        fhit["$channel/jlhit/eventnumber"] = dsp_pars.eventnumber
        fhit["$channel/jlhit/timestamp"] = dsp_pars.timestamp
        fhit["$channel/jlhit/qc"] = dsp_pars.qc
        close(fhit)
        println("hit processing done for $filekey")
    end
    println("DONE! To read the generated hit files: read_ldata(data, :jlhit, category, period, run, channel)")
end
export process_hit