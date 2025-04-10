function rms(x::Vector{T}) where {T <: Real} 
    sqrt(mean(x.^2))
end
export rms
"""
    filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig, τ_pz::Quantity{T}; ft::Quantity{T}= 0.0u"µs", τ_cusp::Quantity{<:AbstractFloat} = 10000000.0u"µs", τ_zac::Quantity{<:AbstractFloat} = 10000000.0u"µs" ) where T<:Real
    filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig; kwargs... )
DSP filter optimization to find best rise-time (for a given flat-top time) to minimize ENC noise. This is an alternative way to calculate the enc noise compared to `dsp_trap_rt_optimization`.
Strategy: 
    filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig, τ_pz::Quantity{T}; ft::Quantity{T}= 0.0u"µs", τ_cusp::Quantity{<:AbstractFloat} = 10000000.0u"µs", τ_zac::Quantity{<:AbstractFloat} = 10000000.0u"µs" ) where T<:Real
- Shift waveforms to have a common baseline, and deconvolute them with the pole-zero correction (in case τ_pz > 0.0u"µs" )
- Filter waveforms with given rise-time and flat-top time
- Build histogram out of all samples in baseline from all waveforms. Remove bins at the beginning and end of the waveform to avoid edge effects.
- Calculate the RMS of the baseline noise --> ENC noise

Inputs:
- `filter_type::Symbol`: filter type (:trap or :cusp)
- `wvfs::ArrayOfRDWaveforms`: raw waveforms to be filtered
- `dsp_config::DSPConfig`: DSP configuration object containing relevant parameters: `ft`, `grid_rt`, `bl_window`, `flt_length_cusp`,  `flt_length_zac`
- `τ_pz::Quantity{T}`: pole-zero decay time. If τ_pz = 0.0u"µs" (or none given), no pole-zero correction is applied.
Optional inputs:
- `ft::Quantity{T}`: fixed flat-top time for optimization
- `τ_cusp::Quantity{<:AbstractFloat}`: cusp decay time; only relevant for cusp filter
- `τ_zac::Quantity{<:AbstractFloat}`: cusp decay time; only relevant for zac filter
"""
function filteropt_rt_optimization_blnoise end
export filteropt_rt_optimization_blnoise
function filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig, τ_pz::Quantity{T}; 
             ft::Quantity{T}= 0.0u"µs", τ_cusp::Quantity{<:AbstractFloat} = 10000000.0u"µs", τ_zac::Quantity{<:AbstractFloat} = 10000000.0u"µs" ) where T<:Real

    # gather config parameters 
    grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
    if ft == 0.0u"µs"
        (_, ft) = get_fltpars(PropDict(), :trap, dsp_config)
    end
    bl_window = dsp_config.bl_window

    # cusp-specific filter parameters 
    flt_length_cusp = getproperty(dsp_config, Symbol("flt_length_cusp"))
    cusp_scale = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # zac-specific filter parameters 
    flt_length_zac = getproperty(dsp_config, Symbol("flt_length_zac"))
    zac_scale = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # waveforms: shift baseline and deconvolute (pole-zero)
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
    wvfs_shift = shift_waveform.(wvfs, -bl_stats.mean)
    wvfs_pz = if τ_pz > 0.0u"µs"
        deconv_flt = InvCRFilter(τ_pz)
        deconv_flt.(wvfs_shift)
    else
        @debug "no decay time given - skip pole-zero correction"
        wvfs_shift
    end 
    # STEP 1: rise-time optimization 
    # calculate baseline rms after filtering 
    # filtered waveforms are already truncated to valid windows. no additional truncation needed.
    noise = zeros(length(grid_rt))
    for (i, rt) in enumerate(collect(grid_rt))
        if 2 * rt + ft >= wvfs_pz[1].time[end]
            noise[i] = NaN
        else
            if filter_type == :trap
                wvfs_flt =  TrapezoidalChargeFilter(rt, ft).(wvfs_pz)   # filtered waveforms 
            elseif filter_type == :cusp
                wvfs_flt = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale).(wvfs_pz) 
            elseif filter_type == :zac
                wvfs_flt = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale).(wvfs_pz)
            end
            valid_bins = findall(leftendpoint(bl_window) .<=  wvfs_flt[1].time .<=rightendpoint(bl_window)) # bins within baseline of filtered waveform 
            bl_trap = filter.(isfinite, map(x-> x.signal[valid_bins], wvfs_flt)) # baselines - cut bins in beginning and end 
            noise[i] = if length(valid_bins) > 0
                rms(vcat(bl_trap...))
            else
                NaN
            end
        end 
    end

    # find optimal shaping time: 
    if  all(.!isfinite.(noise))
        @error "Error: no finite noise values found for filter type $filter_type - try adjusting grid_rt."
    end

    # remove all NaN values from noise array  (too long rise times)
    grid_rt = ustrip.(collect(grid_rt))[findall(isfinite.(noise))]
    noise = noise[findall(isfinite.(noise))]

    result, report = let noise = noise, grid_rt = grid_rt
        f_interp = if length(noise) >= 4
            BSinterpolate(grid_rt, noise, BSplineOrder(4))
        elseif length(noise) >= 2
            LinearInterpolation(grid_rt, noise)
        else 
            x -> NaN
        end
        result = optimize(f_interp, minimum(grid_rt), maximum(grid_rt)) 
        rt_opt = Optim.minimizer(result)*u"µs" # optimial rise time for trap filter
        min_noise = Optim.minimum(result)
        (rt_opt = rt_opt, min_noise = min_noise), (rt_opt = rt_opt, min_noise = min_noise, rt = grid_rt.*u"µs", noise = noise, f_interp = f_interp, ft = ft)
    end
    return result, report 
end

function noise_sweep(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig; kwargs... ) 
    result, report  = filteropt_rt_optimization_blnoise(filter_type, wvfs, dsp_config,  0.0u"µs"; kwargs... )

    return (rt_opt = result.rt_opt, min_noise = result.min_noise, rt = report.rt, noise = report.noise, ft = report.ft), report 
end 
export noise_sweep




