using HDF5
using Distances
using Statistics, LinearAlgebra
using MLJ, MLJClusteringInterface
using CairoMakie, LegendMakie, Makie
using LegendHDF5IO 
using LegendDataManagement
using RadiationDetectorDSP, RadiationDetectorSignals
using LegendDSP
using IntervalSets
using TypedTables: Table as TTable 
using PropDicts
using ColorSchemes
using Printf
using StatsBase 
using Juleanita
using Unitful
using Random
include("$(@__DIR__)/utils_ml.jl")

# data settings and config 
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(52)
channel = ChannelId(1)
category = DataCategory(:cal)
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)

# ml settings: 
pars_ml = try 
    @info "Load ML for $category-$period-$run-$channel"
    asic.par[category].rpars.qc_ml[period,run, channel]
catch 
    PropDict()
end 

if haskey(pars_ml, :ap_opt)
    @info "use optimized AP hyperparameters"
    preference_quantile = pars_ml.ap_opt.opt.preference_quantile
    damp = pars_ml.ap_opt.opt.damp
    maxiter = pars_ml.ap_opt.maxiter
    tol = pars_ml.ap_opt.tol

    # use same random samples as for hyperparameter optimization
    _eventnumber = read_ldata((:eventnumber), asic, DataTier(:raw), filekeys, channel)
    _idx = findall(map(x -> x in pars_ml.ap_opt.wvfs_train_eventnumber, _eventnumber))
else
    @info "hyperparameter not optimized. Use best-guess AP hyperparameters:"
    preference_quantile = 0.01#3
    damp = 0.99
    maxiter = 500
    tol = 1.0e-6

    # draw random samples for the training set
    wvf_max = Int(1e5)
    nsamples = 10000
    rng = MersenneTwister(1234)
    _idx = randperm(rng, wvf_max)[1:nsamples]
end

# load waveforms 
data_raw = TTable(read_ldata(asic, DataTier(:raw), filekeys, channel))[_idx]
wvfs_train_raw = data_raw.waveform
wvfs_train_eventnumber = data_raw.eventnumber

# baseline-shift and normalize waveforms 
wvfs_train = normalize_waveforms(wvfs_train_raw, dsp_config.bl_window)

# sanity plot 5 random waveforms 
# plot settings 
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekeys[1], :ml_qualitycuts) * "/"
for _ in 1:5
    let i = rand(1:length(wvfs_train)) 
        fig = Figure()
        ax = Axis(fig[1, 1], title = "Waveform $i", xlabel =  "Time ($(unit(wvfs_train_raw[i].time[1])))", ylabel = "Signal (a. u.)")
        lines!(ax, ustrip.(wvfs_train_raw[i].time), wvfs_train_raw[i].signal, label = "original")
        lines!(ax, ustrip.(wvfs_train[i].time), wvfs_train[i].signal, label = "normalized")
        axislegend(position = :lt)
        
        _plot_path = joinpath(plt_folder, "waveforms/")
        if !isdir(_plot_path)
            mkpath(_plot_path)
        end
        _pname = "$(_plot_path)AP_waveforms_normalized_$i.png"
        save(_pname, fig)
        @info "Save waveform plot to $(_pname)"
        fig
    end
end 

# run affinity propagation
result_ap, report_ap = trainAP(wvfs_train; 
    preference_quantile = preference_quantile, 
    damp = damp, 
    maxiter = maxiter, 
    tol = tol);
# add waveform eventnumbers to result; 
result_ap = merge(result_ap, (waveforms = merge(result_ap.waveforms, (train_eventnumber = wvfs_train_eventnumber,)),))

if result_ap.ap.ncluster > 120
    @error "Affinity propagation clustering resulted in $(result_ap.ap.ncluster) clusters. This is more than 100 clusters. Modify preference_quantile and damp!"
end

# plot all cluster centers: exemplars with ap labels 
plot_col = get(ColorSchemes.tol_muted, range(0.0, 1.0, length=report_ap.ap.ncluster));
fig_ex_ap = plot_APexemplars(report_ap.exemplars.centers, string.(report_ap.exemplars.labels), plot_col)
_plot_path = joinpath(plt_folder, "AP_exemplars/")
if !isdir(_plot_path)
    mkpath(_plot_path)
end
_pname = "$(_plot_path)/AP_exemplars_damp$(report_ap.ap.damp)_qpref$(report_ap.ap.preference_quantile).png"
save(_pname, fig_ex_ap)
@info "Save exemplars plot to $(_pname)"
fig_ex_ap

# save intermediate results from AP 
pars_ml[:ap] = result_ap
writelprops(asic.par[category].rpars.qc_ml[period], run, PropDict("$channel" => pars_ml))
@info "Save ML-AP results to pars (intermediate,  labelling not done yet)"
