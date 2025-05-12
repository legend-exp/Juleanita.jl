# to do: optimize hyperparameters for the AffinityPropagation
# - preference 
# - damp
using HDF5
using Distances
using Dates
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
using Base.Threads
include("$(@__DIR__)/utils_ml.jl")

date = string(today())

# data settings and config 
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(52)
channel = ChannelId(1)
category = DataCategory(:cal)
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)

# ml settings: 
maxiter = 1000
tol = 1.0e-6
npreference = 10
ndamp = 5
preference_quantiles = vec(hcat(range(0.01, 0.05, 5)...,range(0.1, 0.5, npreference-5)...))
dampings = range(0.5, 0.99, ndamp) # values from Esteban (0.85..0.99). otherwise high probability that AP does not converge. 

# load waveforms and prepare for AP
wvf_max = Int(1e5)
nsamples = 10000
rng = MersenneTwister(1234)
_idx = randperm(rng, wvf_max)[1:nsamples]

data_raw = TTable(read_ldata(asic, DataTier(:raw), filekeys, channel))[_idx]
wvfs_train_raw = data_raw.waveform
wvfs_train_eventnumber = data_raw.eventnumber

# baseline-shift and normalize waveforms 
wvfs_train = normalize_waveforms(wvfs_train_raw, dsp_config.bl_window)

# transform waveforms to a matrix of size (n_waveforms, n_samples) for ML clustering
wvfs_matrix = transpose(hcat(wvfs_train.signal...))
if size(wvfs_matrix) !== (length(wvfs_train), length(wvfs_train[1].signal)) 
    error("Waveform matrix size is not correct - it should be (n_waveforms, n_samples), but is $(size(wvfs_matrix))")
end
X_train = MLJ.table(wvfs_matrix)

# Similary matrix
S =  - pairwise(Cityblock(), wvfs_matrix, dims = 1)
Su = S[triu(trues(size(S)))]
preferences = quantile.(Ref(Su), preference_quantiles)

# get number of cluster over grid
grid_clusters = fill(0, length(preferences), length(dampings))
grid_converged = fill(false, length(preferences), length(dampings))
grid_iterations = fill(0, length(preferences), length(dampings))
_, time_grid, _, _, _ = @timed Threads.@threads for idx in 1:(length(preferences) * length(dampings))
    p = div(idx - 1, length(dampings)) + 1
    d = mod(idx - 1, length(dampings)) + 1
    _damp = dampings[d]
    _preference = preferences[p]

    local model = AffinityPropagation(
            damp = _damp, 
            maxiter = maxiter, 
            tol = tol, 
            preference =  _preference, 
            metric = Cityblock())
    local _machine = machine(model)
    MLJ.predict(_machine, X_train)
    local _report = report(_machine)
    grid_clusters[p, d] = length(_report.cluster_labels)
    grid_converged[p, d] = _report.converged
    grid_iterations[p, d] = _report.iterations
end

# find hyperparameters that give close to 100 clusters and have converged
_clusters = copy(grid_clusters)
_clusters[.!grid_converged] .= 1e8 # ignore not-converged clusters 
min_diff = minimum(abs.(_clusters .-100))
_idx_opt = findall(abs.(_clusters .-100) .== min_diff)
_idx_opt = _idx_opt[argmin(grid_iterations[_idx_opt])] # take the one with less iterations 

opt = (
    preference_quantile = preference_quantiles[_idx_opt[1]],
    preference = preferences[_idx_opt[1]],
    damp = dampings[_idx_opt[2]],
    nclusters = _clusters[_idx_opt],
)

ap_opt = (
    grid_clusters = grid_clusters,
    grid_converged = grid_converged,
    grid_iterations = grid_iterations,
    time_grid = time_grid,
    preferences = preferences,
    preference_quantiles = preference_quantiles,
    dampings = dampings,
    wvfs_train_eventnumber = wvfs_train_eventnumber,
    maxiter = maxiter,
    tol = tol,
    date = date,
    opt = opt
)

# save optimization results from AP 
pars_ml = PropDict(:ap_opt => ap_opt)
writelprops(asic.par[category].rpars.qc_ml[period], run, PropDict("$channel" => pars_ml))

# create plot and save a report. 
pars_ml = asic.par[category].rpars.qc_ml[period, run, channel]
preference_quantiles = pars_ml.ap_opt.preference_quantiles
dampings = pars_ml.ap_opt.dampings
clusters = hcat(pars_ml.ap_opt.grid_clusters...)
converged = hcat(pars_ml.ap_opt.grid_converged...)
clusters[.!converged] .= -1 # remove not-converged clusters from plot 
iterations = hcat(pars_ml.ap_opt.grid_iterations...)
maxiter = pars_ml.ap_opt.maxiter
opt = pars_ml.ap_opt.opt
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekeys[1], :ml_qualitycuts) * "/AP_hyperpars/"
if !isdir(plt_folder)
    mkpath(plt_folder)
end

function plot_matrix(matrix::T; label = "clusters") where T<:AbstractMatrix
    x = 1:length(preference_quantiles)
    fig = Figure()
    ax = Axis(fig[1,1], title ="Affinity propagation -  $(label)",
                        xlabel = "Preference quantile", ylabel = "Damping",  
                        xticks = (x, string.(round.(preference_quantiles, digits=2))),
                        yticks = (dampings, string.(round.(dampings, digits=2))),
                        titlesize = 14)
    if any(clusters.< 0)
        col_range = (0, maximum(matrix))
    else
        col_range = (minimum(matrix), maximum(matrix))
    end                  
    heatmap!(ax, x, dampings, matrix, colormap = :isoluminant_cm_70_c39_n256, lowclip = :lightgrey, colorrange = col_range)
    # Overlay text annotations for each cell
    for i in eachindex(preference_quantiles)
        for j in eachindex(dampings)
            if matrix[i, j] == -1
                continue # skip not-converged cells
            end
            text!(
                ax,
                x[i],  # x-coordinate (preference value)
                dampings[j],        # y-coordinate (damp value)
                text = string(matrix[i, j]),  # Value in the cell
                align = (:center, :center),     # Center the text in the cell
                fontsize = 12,                  # Adjust font size as needed
                color = :white                  # Text color
            )
        end
    end
    fig
    # Save the figure
    pname = "$(plt_folder)AP_hyperpars_opt_$(replace(label, " " => ""))_$(length(dampings))dampings$(minimum(dampings))-$(maximum(dampings))_$(length(preference_quantiles))qprefs$(minimum(preference_quantiles))-$(maximum(preference_quantiles)).png"
    save(pname, fig)
    @info "Figure saved to $pname"
    fig
end
fig_iterations = plot_matrix(iterations; label = "Number of iterations (max $(maxiter))")
fig_cluster = plot_matrix(clusters; label = "Number of clusters")
fig_converged = plot_matrix(converged; label = "Convergence")

# report later saved to this folder
rname = "$(plt_folder)AP_hyperpars_opt_$(length(dampings))dampings$(minimum(dampings))-$(maximum(dampings))_$(length(preference_quantiles))qprefs$(minimum(preference_quantiles))-$(maximum(preference_quantiles)).md"   
# report 
begin 
    _lreport = lreport()
    lreport!(_lreport, "# Affinity Propagation hyperparameters optimization")
    lreport!(_lreport, "- $date")
    lreport!(_lreport, "- $category-$period-$run-$channel")
    lreport!(_lreport, "- Computing time : $(round(pars_ml.ap_opt.time_grid/3600, digits = 2)) hours")
    lreport!(_lreport, "We perfrom AP  over a grid of different dampings and preference quantiles. ")
    lreport!(_lreport, "The Training is based on $(length(pars_ml.ap_opt.wvfs_train_eventnumber)) randomly drawn waveforms (no repetitions)")
    lreport!(_lreport, "## Optimal hyperparameters:")
    lreport!(_lreport, "Found *optimal* hyperparameter combination for which AP:")
    lreport!(_lreport, "- converges")
    lreport!(_lreport, "- gives closest to 100 clusters")
    lreport!(_lreport, "- has smallest number of iterations (in case serveral APs with same minimal number of clusters)")
    lreport!(_lreport, "- **result: damping = $(opt.damp), preference = $(opt.preference_quantile) (quantile) $(round(opt.preference, digits = 1)) (abs.), number of clusters = $(opt.nclusters)**" )
    lreport!(_lreport, "## 1. Number of clusters/exemplars:")
    lreport!(_lreport, "The following figure shows the number of resulting cluster (exemplars) for different combinations of damping and preference. The optimal (damp, preference)-combination should have about 100 exemplars. ")
    lreport!(_lreport, "Grid points for which AP did not converge are marked as grey.")
    lreport!(_lreport, fig_cluster)
    lreport!(_lreport, ". ")
    lreport!(_lreport, "## 2. Convergence:")
    lreport!(_lreport, fig_converged )
    lreport!(_lreport, ". ")
    lreport!(_lreport, "## 3. Number of iterations:")
    lreport!(_lreport, "The maximum number of iterations is set to $(maxiter). After that, AP stops whether converged or not.")
    lreport!(_lreport, fig_iterations)
    # column_names = vcat([Symbol("pref$(round(p, digits=2))") for p in preference_quantiles]...)
    # tab_clusters = TTable(damping = dampings; (column_names .=> Vector(eachrow(clusters)))...)
    # tab_converged = TTable(damping = dampings; (column_names .=> Vector(eachrow(converged)))...)
    # lreport!(_lreport, tab_clusters)
    # lreport!(_lreport, tab_converged)
    writelreport(rname, _lreport)
    @info "Report saved to $rname"
    _lreport
end 

