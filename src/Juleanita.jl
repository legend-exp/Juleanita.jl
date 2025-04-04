module Juleanita

using BSplineKit: BSplineOrder, interpolate as BSinterpolate
using CSV
using Dates 
using DataFrames
using HDF5
using IntervalSets
using JSON
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using LegendDSP
using LegendDSP: get_fltpars
using LegendSpecFits
using LegendSpecFits: get_friedman_diaconis_bin_width
using LinearAlgebra
using Makie, CairoMakie, LegendMakie
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
export mvalue, muncert
using Measures
using Optim
using Printf
using PropDicts 
using RadiationDetectorDSP
using RadiationDetectorSignals
using RadiationSpectra
using StatsBase
using TypedTables
using Unitful

include("apply_qc.jl")
include("filteropt_rt_optimization_blnoise.jl")
include("fit_linearity.jl")
include("IO_csv.jl")
include("utils.jl")
include("utils_physics.jl")
include("simple_dsp.jl")
include("process_peak_split.jl")
include("process_decaytime.jl")
include("process_filteropt.jl")
include("process_ctc.jl")
include("process_dsp.jl")
include("process_qualitycuts.jl")
include("process_energy_calibration.jl")
include("process_hit.jl")
include("process_peakfits.jl")
end

