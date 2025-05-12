using LegendDataManagement
using RadiationDetectorSignals, RadiationDetectorDSP
using Unitful
using IntervalSets
using PropDicts
using StatsBase
using Makie, CairoMakie, LegendMakie
"""
    normalize_waveforms(wvfs::ArrayOfRDWaveforms, bl_window::ClosedInterval{<:Quantity})
Normalize the input waveforms by baseline-shifting and normalizing to the abs. maximum amplitude.
"""
function normalize_waveforms(wvfs::ArrayOfRDWaveforms, bl_window::ClosedInterval{<:Quantity})
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
    wvfs_shift = shift_waveform.(wvfs, -bl_stats.mean)
    wvf_absmax = maximum.(map(x-> abs.(x), wvfs_shift.signal))
    return multiply_waveform.(wvfs_shift, 1 ./ wvf_absmax)
end

"""
    trainAP(wvfs_train::ArrayOfRDWaveforms; preference_quantile::T = 0.5, damp::T = 0.1, maxiter::Int = 200, tol::T = 1.0e-6) where T<:Real
Train the AffinityPropagation model on the input waveforms.
"""
function trainAP(wvfs_train::ArrayOfRDWaveforms; preference_quantile::T = 0.5, damp::T = 0.1, maxiter::Int = 200, tol::T = 1.0e-6) where T<:Real
    # transform waveforms to a matrix of size (n_waveforms, n_samples) for ML clustering
    wvfs_matrix = transpose(hcat(wvfs_train.signal...))
    if size(wvfs_matrix) !== (length(wvfs_train), length(wvfs_train[1].signal)) 
        error("Waveform matrix size is not correct - it should be (n_waveforms, n_samples), but is $(size(wvfs_matrix))")
    end
    X = MLJ.table(wvfs_matrix)

    # Similary matrix
    S =  - pairwise(Cityblock(), wvfs_matrix, dims = 1)
    Su = S[triu(trues(size(S)))]
    preference = quantile(Su, preference_quantile)

    # prepare model
    model = AffinityPropagation(
                damp = damp, 
                maxiter = maxiter, 
                tol = tol, 
                preference = preference, 
                metric = Cityblock())

    # apply model to data
    begin
        _machine = machine(model)
        waveform_ap_labels = collect(MLJ.predict(_machine, X)) # these are the labels of each waveform (assigned to a cluster)
    end        

    # evaluate the model 
    _report = MLJ.report(_machine)

    # check convergence 
    if !_report.converged
        @warn "AffinityPropagation did not converge. Check the model parameters."
    end
    @info "$(length(_report.cluster_labels)) clusters found"

    # summarize results
    result_ap = (
        waveforms = (ap_labels = waveform_ap_labels, ),
        ap = (damp = damp,  preference = preference, preference_quantile = preference_quantile, 
            iterations = _report.iterations, ncluster = length(_report.exemplars), converged = _report.converged),
        exemplars = (idx = _report.exemplars, eventnumber =  wvfs_train_eventnumber[_report.exemplars],
                    labels = _report.cluster_labels),
    )

    report_ap = (waveforms = result_ap.waveforms,
                    ap = result_ap.ap,
                    exemplars = merge(result_ap.exemplars, (centers = wvfs_train[_report.exemplars],) ))

    return result_ap, report_ap 
end 

"""
    plot_APexemplars(centers::ArrayOfRDWaveforms, labels::Vector, colors::Vector)
Plot the cluster centers (exemplars) of the AffinityPropagation clustering.

"""
function plot_APexemplars(centers, labels::Vector, colors::Vector)
    ncluster = length(centers)
    ncol = round(Int, sqrt(ncluster))
    nrow = ceil(Int, ncluster/ncol)
    fig = Figure(size=(ncol*100,nrow*110), figure_padding = 5)
	for i in 1:ncluster
        # wf_idx = result_ap.exemplars.idx[i]
        # qc_label = result_ap.waveforms.qc_labels[wf_idx]
        _row = div(i-1,ncol)+1
        _col = mod(i-1,ncol)+1
        _ax =  Axis(fig[_row,_col], xticklabelsize = 5, yticklabelsize = 5)
		lines!(_ax, ustrip.(centers[i].time), centers[i].signal, linewidth = 1.5, color = colors[i])
        Label(fig[_row, _col], "$(labels[i])", fontsize = 10, padding = (10, 10),  tellwidth = false)
        hidedecorations!(_ax)
	end
    Label(fig[0, :], "Cluster Center Exemplar Waveforms", fontsize = 20, tellwidth = false)
	rowgap!(fig.layout, 5)
    colgap!(fig.layout, 5)
    return fig
end

function ml_filename(data::LegendData, category::DataCategory, period::DataPeriod, run::DataRun )
    ml_folder = data.tier[DataTier(:jlqcml), category , period, run] * "/"
    if !ispath(ml_folder)
        mkpath(ml_folder)
    end
    filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
    ml_folder * string(filekey) * "-tier_jlqcml.jld2"
end

"""
    dwt(waveforms::ArrayOfRDWaveforms; nlevels::Int = 2)
Discrete wavelet transform (DWT) of the input waveforms using Haar filter
"""
function dwt(waveforms::ArrayOfRDWaveforms; nlevels::Int = 2)
    # create Haar filter 
    haar_flt = HaarAveragingFilter(2)  

    # Apply Haar filter nlevels times
    wvfs_flt = copy(waveforms)
    for _ in 1:nlevels
        wvfs_flt = haar_flt.(wvfs_flt)
    end

    # normalize with max of absolute extrema values
    norm_fact = map(x -> max(abs(first(x)), abs(last(x))), extrema.(wvfs_flt.signal))
    replace!(norm_fact, 0.0 => one(first(norm_fact)))
    wvfs_flt = multiply_waveform.(wvfs_flt, 1 ./ norm_fact)
end 

"""
    plot_SVM_QCeff(label_pred::Vector, pars_mlap::PropDict; dataset = "", accuracy = NaN)
Plot the distribution of qc-labels based on the trained full AP-SVM model
If accuracy is provided, it will be shown in the plot legend 
"""
function plot_SVM_QCeff(label_pred::Vector, pars_mlap::PropDict; dataset = "", accuracy = NaN, plot_name = NaN)
    _qc_cat =  parse.(Int, string.(keys(countmap(label_pred))))
    qc_cat = [pars_mlap.legend.qc_labels[_qc_cat[i]] for i in eachindex(_qc_cat)]
    x  = collect(1:length(qc_cat))
    y  = parse.(Int, string.(values(countmap(label_pred)))) ./ length(label_pred)
    label_plot = if !isfinite(accuracy)
       "$(length(label_pred)) waveforms ($dataset)"
    else 
       "$(length(label_pred)) waveforms ($dataset) \naccuracy: $(round(accuracy, digits = 3))"
    end

    fig = Figure()
    ax = Axis(fig[1, 1], 
        limits = ((nothing, nothing), (0, 1.05)),
        xlabel = "QC category", 
        ylabel = "Fraction of training waveforms", 
        title = title = Juleanita.get_plottitle(filekeys[1], _channel2detector(asic, channel), "AP-SVM Quality cuts"),
        titlesize = 16)
    ax.xticks = x
    ax.xtickformat = x -> qc_cat
    barplot!(ax, x, y, bar_labels = :y, label = label_plot)
    axislegend()
    if isa(plot_name, String)
        save(plot_name , fig)
        @info "Save SVM efficiency plot to $(plot_name)"
    end 
    fig 
end 

