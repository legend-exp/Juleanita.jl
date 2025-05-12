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
include("$(@__DIR__)/utils_ml.jl")

# data settings  
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(52)
channel = ChannelId(1)
category = DataCategory(:cal)

# plot settings
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekeys[1], :ml_qualitycuts) * "/"

# load AP training results (step 01)
ml_file = joinpath(data_path(asic.par[category].rpars.qc_ml[period]), "$run.json" )
pars_ml  = if isfile(ml_file) 
    @info "Load AP training results for $category-$period-$run-$channel"
    asic.par[category].rpars.qc_ml[period,run, channel]
else 
    error("AP training needs to happen before relabelling! Cannot proceed")
end 
result_ap = pars_ml.ap
if haskey(result_ap, :qc_labels) 
    @info "Relabelling already done. If you want to redo relabelling, run following script"
end 

# load waveforms used for AP training
eventnumber = read_ldata((:eventnumber), asic, DataTier(:raw), filekeys, channel)
data_raw = TTable(read_ldata(asic, DataTier(:raw), filekeys, channel))[findall(map(x -> x in result_ap.waveforms.train_eventnumber, eventnumber))]
wvfs_train_raw = data_raw.waveform
@assert data_raw.eventnumber == result_ap.waveforms.train_eventnumber

# plot all cluster centers: exemplars with ap labels 
centers = wvfs_train_raw[result_ap.exemplars.idx]
plot_col = get(ColorSchemes.tol_muted, range(0.0, 1.0, length=result_ap.ap.ncluster));
fig_ex_ap = plot_APexemplars(centers, string.(Int.(result_ap.exemplars.labels)), plot_col)
fig_ex_ap

# THIS STEP HAS TO BE DONE MANUALLY!!!!!
# Rename exemplars/clusters AP-labels ---> LEGEND QC-labels. 1 QC-label can be assigned to more than one AP-label
# re-label the all waveforms by hand. several cluster can belong to the the same QC label. 
@info "Renaming clusters into LEGEND qc-labels. THIS HAS TO BE DONE MANUALLY!!!! "
@info " here is the legend for thr QC labels:"
qc_labels = dataprod_config(asic).qc_ml(filekeys[1]).qc_labels
relabel_nt = (normal = [4,8,9,10,11,12,14, 18, 21, 22, 24, 25, 26, 27, 29, 30,34,35,40,41,42, 44, 46,48, 50,51,52,53, 54, 55, 56,62, 64,66,67,71,72,73,74,75,76,77,78,85, 87, 90, 92, 94,95,97,98,99,102, 103, 104],
    neg_go = [],
    up_slo = [],
    down_slo = [],
    spike = [],
    x_talk = [],
    slo_ri = [],
    early_tr = [],
    late_tr = [],
    sat = [],
    soft_pi = [16,20,63],
    hard_pi = [1,5 , 6, 13, 15,17,19, 23, 28, 31, 32,33,36,37,38,43,45,49,57,58, 59,60, 65, 68, 69,70,79,80,81,82,83, 86,88,89, 93,96,100,101],
    bump = [2, 3, 7, 39, 47, 61, 84, 91], # different than in Legend: here bump in rising edge
    noise = [],
)

# sanity check for relabeling 
begin
    @assert length(vcat(values(relabel_nt)...)) == result_ap.ap.ncluster  "Relabelling does not match number of clusters"
    @assert length(unique(vcat(values(relabel_nt)...))) == length(vcat(values(relabel_nt)...)) "Relabelling has duplicates" 
    if length(unique(vcat(values(relabel_nt)...))) !== result_ap.ap.ncluster
        findall(map(x -> x ∉ vcat(values(relabel_nt)...), range(result_ap.ap.ncluster) ) )
    end
end 

"""
    assign_qc_labels(ap_labels, relabel_nt::NamedTuple)
Convert Vector of AP labels (number between 0 and ncluster) to a Vector of QC labels (0-13).

INPUT: 
- ap_labels::Vector of AP-labels (number between 0 and ncluster)
- relabel_nt::NamedTuple that maps which AP-labels belong to qhich QC-label (string format)
OUTPUT:
- qc_labels: Vector of QC-labels (0-13)
"""
function assign_qc_labels(ap_labels, relabel_nt::NamedTuple)
    qc_labels = fill(0, length(ap_labels))
    for (i, ap_label) in enumerate(ap_labels)
        if ap_label in relabel_nt.normal 
            qc_labels[i] = 0
        elseif ap_label in relabel_nt.neg_go
            qc_labels[i] = 1
        elseif ap_label in relabel_nt.up_slo
            qc_labels[i] = 2
        elseif ap_label in relabel_nt.down_slo
            qc_labels[i] = 3
        elseif ap_label in relabel_nt.spike
            qc_labels[i] = 4
        elseif ap_label in relabel_nt.x_talk
            qc_labels[i] = 5
        elseif ap_label in relabel_nt.slo_ri
            qc_labels[i] = 6
        elseif ap_label in relabel_nt.early_tr
            qc_labels[i] = 7
        elseif ap_label in relabel_nt.late_tr
            qc_labels[i] = 8
        elseif ap_label in relabel_nt.sat
            qc_labels[i] = 9
        elseif ap_label in relabel_nt.soft_pi
            qc_labels[i] = 10
        elseif ap_label in relabel_nt.hard_pi
            qc_labels[i] = 11
        elseif ap_label in relabel_nt.bump
            qc_labels[i] = 12
        elseif ap_label in relabel_nt.noise
            qc_labels[i] = 13
        else
            error("Unknown AP label: $ap_label")
        end
    end 
    return qc_labels
end

# apply relabel-logic to each waveform and exemplars
waveform_qc_labels = assign_qc_labels(result_ap.waveforms.ap_labels, relabel_nt)
exemplar_qc_labels = assign_qc_labels(result_ap.exemplars.labels, relabel_nt)

# add qc_labels in result_ap and save results to pars
waveforms = merge(result_ap.waveforms, PropDict(:qc_labels => waveform_qc_labels))
exemplars = merge(result_ap.exemplars, PropDict(:qc_labels => exemplar_qc_labels))
result_ap = PropDict(:ap => result_ap.ap,
                    :waveforms => merge(result_ap.waveforms, PropDict(:qc_labels => waveform_qc_labels)), 
                    :exemplars => merge(result_ap.exemplars, PropDict(:qc_labels => exemplar_qc_labels)),
                    :legend => PropDict(:qc_labels => qc_labels, :relabel_nt => relabel_nt))
pars_ml[:ap] = result_ap
writelprops(asic.par[category].rpars.qc_ml[period], run, PropDict("$channel" => pars_ml))
@info "Save ML-AP relabel results to pars"


# SANITY PLOT: plot all cluster centers: exemplars with new qc labels
begin 
    plot_qclabel = result_ap.exemplars.qc_labels
    plot_qclabel_str = map(x -> "$(result_ap.legend.qc_labels[x]) ($x)", plot_qclabel) 
    colors = get(ColorSchemes.tol_muted, range(0.0, 1.0, length=14))
    plot_col = colors[plot_qclabel .+ 1 ]
    fig_ex_qc = plot_APexemplars(centers, plot_qclabel_str, plot_col)
    _plot_path = joinpath(plt_folder, "AP_exemplars/")
    if !isdir(_plot_path)
        mkpath(_plot_path)
    end
    _pname = "$(_plot_path)/AP_exemplars_QClabels_damp$(result_ap.ap.damp)_qpref$(result_ap.ap.preference_quantile).png"
    save(_pname, fig_ex_qc)
    @info "Save exemplars with qc labels plot to $(_pname)"
    fig_ex_qc
end


# EFFICIENCY plot figure 
ap_pars = asic.par[category].rpars.qc_ml[period, run, channel].ap
_qc_cat = collect(keys(countmap(ap_pars.waveforms.qc_labels)))
qc_cat = [ap_pars.legend.qc_labels[_qc_cat[i]] for i in eachindex(_qc_cat)]
x  = collect(1:length(qc_cat))
y  = parse.(Int, string.(values(countmap(ap_pars.waveforms.qc_labels)))) ./ length(ap_pars.waveforms.qc_labels)

fig = Figure()
ax = Axis(fig[1, 1], 
    limits = ((nothing, nothing), (0, 1)),
    xlabel = "QC category", 
    ylabel = "Fraction of training waveforms", 
    title = title = Juleanita.get_plottitle(filekeys[1], _channel2detector(asic, channel), "Affinity Propagation Clustering"))
ax.xticks = x
ax.xtickformat = x -> qc_cat
barplot!(ax, x, y, bar_labels = :y, label = "$(length(ap_pars.waveforms.ap_labels)) waveforms")
axislegend()
fig
_plot_path = joinpath(plt_folder, "efficiency/")
if !isdir(_plot_path)
    mkpath(_plot_path)
end
_pname = "$(_plot_path)AP_fractions_damp$(ap_pars.ap.damp)_qpref$(ap_pars.ap.preference_quantile).png"
save(_pname , fig)
fig
