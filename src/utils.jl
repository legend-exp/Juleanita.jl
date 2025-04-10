function _quantile_truncfit(x::Vector; qmin::T = 0.02, qmax::T = 0.98)  where T<:Real
    filter!(isfinite, x)
    quantile(skipmissing(x), qmin),  quantile(skipmissing(x), qmax)
end

function _quantile_truncfit(x::Matrix; qmin::T= 0.02, qmax::T = 0.98) where T<:Real
    nrows = size(x)[1]
    xmins = zeros(nrows,1)
    xmaxs = zeros(nrows,1)

    for i = 1:nrows
        if  isempty(filter(.!ismissing, x[i,:])) || all(isnan.(x[i,:]))
            xmins[i], xmaxs[i] = NaN, NaN
        else
            xmins[i], xmaxs[i] = _quantile_truncfit(x[i,:]; qmin = qmin, qmax = qmax) 
        end
    end

   minimum(filter(isfinite, xmins)), maximum(filter(isfinite, xmaxs))
end

function _channel2detector(data::LegendData, channel::ChannelId)
    detectors = [DetectorId(Symbol(data.metadata.hardware.detectors.germanium[k].name)) for k in keys(data.metadata.hardware.detectors.germanium)]
    channels   = [ChannelId(data.metadata.hardware.detectors.germanium[k].channel) for k in keys(data.metadata.hardware.detectors.germanium)]
    detectors[findfirst(x -> x== channel, channels)]
end
export _channel2detector

function _detector2channel(data::LegendData, detector::DetectorId)
    detectors = [DetectorId(Symbol(data.metadata.hardware.detectors.germanium[k].name)) for k in keys(data.metadata.hardware.detectors.germanium)]
    channels = [ChannelId(data.metadata.hardware.detectors.germanium[k].channel) for k in keys(data.metadata.hardware.detectors.germanium)]
    channels[findfirst(x -> x== detector, detectors)]
end
export  _detector2channel

function _get_pltfilename(data::LegendData, filekey::FileKey, channel::ChannelIdLike, process::Symbol)
    pname = split(LegendDataManagement.LDMUtils.get_pltfilename(data, filekey, channel, process),"/")[end]
    d = LegendDataManagement.LDMUtils.get_pltfolder(data, filekey, process)
    ifelse(isempty(readdir(d)), rm(d), nothing )
    return pname 
end
export _get_pltfilename




