using HDF5, JLD2
using Distances
using Statistics, LinearAlgebra
using MLJ, MLJClusteringInterface
using LIBSVM
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
using Unitful 
using Juleanita
include("$(@__DIR__)/utils_ml.jl")

# data settings
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(52)
channel = ChannelId(1)
category = DataCategory(:cal)
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)

# ml-SVM settings: 
dwt_haar_levels = 6 
train_to_test_split = 0.8
svm_cost = 1.0
svm_gamma = 0.5

# 
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekeys[1], :ml_qualitycuts) * "/"

# load AffinityPropagation results (if not available, give error message)
pars_ml = try 
    @info "Load ML for $category-$period-$run-$channel"
    asic.par[category].rpars.qc_ml[period,run, channel]
catch 
    error("Run ML-AP clustering + relabelling for $category-$period-$run-$channel do not exist. Run 01_AP_hperpars_opt.jl, 02_AP_train.jl and 03_AP_relabel.jl first")
end 
if !haskey(pars_ml.ap.exemplars, :qc_labels)
    error("Run ML-AP relabelling for $category-$period-$run-$channel does not exist. Run  02_AP_relabel.jl first")
end 

# load training waveforms that were also used during AP training and normalize them
eventnumber = read_ldata((:eventnumber), asic, DataTier(:raw), filekeys, channel)
data_raw = TTable(read_ldata(asic, DataTier(:raw), filekeys, channel))[findall(map(x -> x in pars_ml.ap.waveforms.train_eventnumber, eventnumber))]
wvfs_train = normalize_waveforms(data_raw.waveform, dsp_config.bl_window)

# do Discrete Wavelet Transform (DWT) on the waveforms
"""
    dwt(waveforms::ArrayOfRDWaveforms; nlevels::Int = 2)
Discrete wavelet transform (DWT) of the input waveforms using Haar filter
"""
function dwt(waveforms::ArrayOfRDWaveforms; nlevels::Int = 2)
    # create Haar filter 
    haar_flt = HaarAveragingFilter(2)  

    # wvfs_flt = waveforms .|> haar_flt .|> haar_flt
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

wvfs_dwt = dwt(wvfs_train; nlevels = dwt_haar_levels)
let fig = Figure()
    i = rand(1:length(wvfs_dwt))
    ax = Axis(fig[1, 1], title = "Discrete wavelet transform ($i): \ndownsampled + noise reduced", xlabel = "Time (µs)", ylabel = "Norm. Signal", limits = ((nothing, nothing), (nothing, 1.35)))
    lines!(ax, ustrip.(wvfs_train[i].time), wvfs_train[i].signal, label = "Original ($(round(Int, length(wvfs_train[i].time)/1e3))k samples)")
    lines!(ax, ustrip.(wvfs_dwt[i].time), wvfs_dwt[i].signal, label = "DWT level $(dwt_haar_levels) ($(round(Int, length(wvfs_dwt[i].time)/1e3))k samples)", linestyle = :dash)
    axislegend(position = :lt)
    save("$(plt_folder)/waveforms/Waveform_DWT_$i.png", fig)
	fig
end
# sanity check: plot 5 random waveforms
# divide the waveforms into training and test sets 
ntrain = round(Int, length(wvfs_dwt)*train_to_test_split)
labels_train = pars_ml.ap.waveforms.qc_labels[1:ntrain]
labels_test = pars_ml.ap.waveforms.qc_labels[ntrain+1:end]
X_train = hcat(wvfs_dwt.signal...)[:, 1:ntrain]
X_test = hcat(wvfs_dwt.signal...)[:, ntrain+1:end]
@assert size(X_train, 2) == length(labels_train) "Number of samples in X must match the length of label"

# Train an SVM model with a radial basis function (RBF) kernel
model = svmtrain(X_train, labels_train; kernel=LIBSVM.Kernel.RadialBasis, cost = svm_cost, gamma = svm_gamma)
model_predict = Base.Fix1(svmpredict, model)

# Predict labels for the training data
pred_labels_train = svmpredict(model, X_train)
accuracy_train = sum(pred_labels_train[1] .== labels_train) / length(labels_train)
println("Accuracy train: $accuracy_train")

# Predict labels for the test data
pred_labels_test = svmpredict(model, X_test)
accuracy_test = sum(pred_labels_test[1] .== labels_test) / length(labels_test)
println("Accuracy test: $accuracy_test")

# save pars (without model, not possible to json )
result_svm = PropDict(:svm => (cost = svm_cost, gamma = svm_gamma, kernel = Int(model.kernel), dwt_haar_levels = dwt_haar_levels), 
                     :train => (train_to_test_split = train_to_test_split, accuracy_train = accuracy_train, labels_pred = pred_labels_train, labels_true = labels_train),
                     :test => (accuracy_test = accuracy_test, labels_pred = pred_labels_test, labels_true = labels_test))
pars_ml = merge(pars_ml, PropDict(:svm => result_svm))
writelprops(asic.par[category].rpars.qc_ml[period], run, PropDict("$channel" => pars_ml))
@info "Save ML-AP-SVM relabel results to pars"
# save everything including trained model to JLD2 file.  
ml_file = ml_filename(asic, category, period, run)
pars_ml = merge(pars_ml, PropDict(:svm => merge(result_svm, PropDict(:model => (model = model, model_pred = model_predict)))))
save(ml_file, Dict("$channel" => pars_ml))

# sanity plots: distribution of predicted labels for training and test data 
_plot_path = joinpath(plt_folder, "SVM/")
if !isdir(_plot_path)
    mkpath(_plot_path)
end
pname = (dataset) -> "$(_plot_path)SVM_$(dataset)_cost$(pars_ml.svm.svm.cost)_gamma$(pars_ml.svm.svm.gamma).png"    
plot_SVM_QCeff(pred_labels_train[1], pars_ml.ap; dataset = "train", accuracy = accuracy_train, plot_name = pname("train"))
plot_SVM_QCeff(pred_labels_test[1], pars_ml.ap; dataset = "test",  accuracy = accuracy_test, plot_name = pname("test"))

