using Juleanita
using PropDicts
using JLD2
using LegendDataManagement
using RadiationDetectorSignals, RadiationDetectorDSP
using LegendDSP
using Unitful
using IntervalSets
using Makie, CairoMakie, LegendMakie
include("$(@__DIR__)/utils_ml.jl")

# data settings and configs
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(52)
channel = ChannelId(1)
category = DataCategory(:cal)
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)

# load ML quality cut trained model from file 
ml_file = ml_filename(asic, category, period, run)
qc_pd = JLD2.load(ml_file)["$channel"]
qc_model_func = qc_pd.svm.model.model_pred 

# load waveforms and prepare for ML model 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
wvfs_raw = read_ldata((:waveform), asic, DataTier(:raw), filekeys, channel)

# baseline-shift and normalize waveforms 
wvfs_train = normalize_waveforms(wvfs_raw, dsp_config.bl_window)
wvfs_dwt = dwt(wvfs_train; nlevels = dwt_haar_levels)

# apply ML model to waveforms
label_pred = qc_model_func(hcat(wvfs_dwt.signal...))
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekeys[1], :ml_qualitycuts) * "/SVM/"

pname = (dataset) -> "$(plt_folder)SVM_test_cost$(pars_ml.svm.svm.cost)_gamma$(pars_ml.svm.svm.gamma).png"    
plot_SVM_QCeff(label_pred[1], qc_pd.ap; dataset = "new", accuracy = NaN, plot_name = pname("new"))
