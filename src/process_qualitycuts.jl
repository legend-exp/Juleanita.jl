"""
    process_qualitycuts(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; reprocess::Bool = false, qc_config::PropDict = data.metadata.config.qc.qc_config.default)
apply quality cuts based on dsp parameters 
inputs: 
    data: LegendData object
    period: DataPeriod object
    run: DataRun object
    category: Symbol or DataCategory object
    channel: ChannelId object
    reprocess: Bool, default false
    qc_config: PropDict, default data.metadata.config.qc.qc_config.default
"""
function process_qualitycuts end
export process_qualitycuts
function process_qualitycuts(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; 
                            reprocess::Bool = false, qc_config::PropDict = data.metadata.config.qc.qc_config.default)

    # check if quality cut pars already exist
    qc_file = joinpath(mkpath(data_path(data.par[category].rpars.qc[period])), "$(string(run)).json")
     if isfile(qc_file) && !reprocess
         @info "Quality cuts (qc) file already exist for $category period $period - run $run - channel $channel - you're done!"
         return
     end

    # load dsp parameters 
    dsp_par = Table(read_ldata(data, :jldsp, category, period, run, channel))

    # calculate quality cuts. defined in src/apply_qc.jl
    qc = apply_qc(dsp_par, qc_config)

    # add event number and timestamp
    qc[:timestamp] = dsp_par.timestamp
    qc[:eventnumber] = dsp_par.eventnumber

    # save results to pars 
    result_qc = PropDict(Dict("$channel" => qc))
    writelprops(data.par[category].rpars.qc[period], run, result_qc)
    @debug "Saved qc pars to disk"

    # add qcflag to dsp tier. 
    filekeys = search_disk(FileKey, data.tier[:jldsp, category , period, run])
    dsp_folder =  data.tier[DataTier(:jldsp), category , period, run] * "/"
    for fk in filekeys
        fname = dsp_folder .* string(fk) .* "-tier_jldsp.lh5"
        h5open(fname, "r+") do f
            if haskey(f, "$channel/jldsp/qc")
                delete_object(f, "$channel/jldsp/qc")
                @debug "remove old qc values in $fname"
            end
        end
        eventnumbers_fk = read_ldata(:eventnumber, data, :jldsp, fk, channel);
        qc_indices = findall(x -> x in eventnumbers_fk, qc.eventnumber)
        fdsp = lh5open(fname, "cw")
        fdsp["$channel/jldsp/qc"] = qc.wvf_keep.all[qc_indices]
        close(fdsp)
       @debug "add qc to dsp file for $fk"
    end 
    return qc 
end