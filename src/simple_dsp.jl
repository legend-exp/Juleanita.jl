"""
    simple_dsp(data::Q, dsp_config::DSPConfig; τ_pz::Quantity{T} = 0.0u"µs", pars_filter::PropDict) where {Q <: Table, T<:Real}
dsp routine which calculates all relevant parameters for the waveform analysis
- `data::Q`: input data, e.g. raw file, peakfile 
- `dsp_config::DSPConfig`: DSP configuration object
- `τ_pz::Quantity{T}`: pole-zero decay time. If τ_pz = 0.0u"µs" (or none given), no pole-zero correction is applied.
- `pars_filter::PropDict`: filter parameters for the different filter types. Use PropDict() to get default values from config.
"""
function simple_dsp(data::Q, dsp_config::DSPConfig; τ_pz::Quantity{T} = 0.0u"µs", pars_filter::PropDict) where {Q <: Table, T<:Real}
    wvfs = data.waveform

    #### get filter parameters from config. if otimized values are not found; use default.  
    trap_rt, trap_ft = mvalue.(get_fltpars(pars_filter, :trap, dsp_config))
    cusp_rt, cusp_ft = mvalue.(get_fltpars(pars_filter, :cusp, dsp_config))
    zac_rt, zac_ft   = mvalue.(get_fltpars(pars_filter, :zac, dsp_config))
    sg_wl            = mvalue(get_fltpars(pars_filter, :sg, dsp_config))

    # get CUSP and ZAC filter length and flt scale
    flt_length_zac              = getproperty(dsp_config, Symbol("flt_length_zac")) 
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))
    flt_length_cusp             = getproperty(dsp_config, Symbol("flt_length_cusp")) 
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"
    τ_zac = 10000000.0u"µs"

    ################## ACTUAL WAVEFORM FILTERING AND RECONSTRUCTION, ANALYSIS BEGINS HERE ##################
    bl_stats = signalstats.(wvfs, leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
   
    # tail analysis 
    tail_stats = tailstats.(wvfs, leftendpoint(dsp_config.tail_window), rightendpoint(dsp_config.tail_window))

    # get raw wvf maximum/minimum
    wvf_max = maximum.(wvfs.signal)
    wvf_min = minimum.(wvfs.signal)
    
    # deconvolute waveform: pole-zero correction. Use pre-defined tau from decay time analysis OR median decay time for all waveforms
    if τ_pz == 0.0u"µs"
        τ_pz = median(filter!(isfinite, tail_stats.τ))
    end
    deconv_flt = InvCRFilter(τ_pz)
    wvfs = deconv_flt.(wvfs)
  
    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs, leftendpoint(dsp_config.tail_window), rightendpoint(dsp_config.tail_window))

    # t0 determination
    t0 = get_t0(wvfs, dsp_config.t0_threshold; flt_pars=dsp_config.kwargs_pars.t0_flt_pars, mintot=dsp_config.kwargs_pars.t0_mintot)
    
    # if all waveforms are saturated set threshold to 1.0 to avoid numerical problems
    # replace!(wvf_max, zero(wvf_max[1]) => one(wvf_max[1]))
    
    # get threshold points in rise
    t10 = get_threshold(wvfs, wvf_max .* 0.1; mintot=dsp_config.kwargs_pars.tx_mintot)
    t20 = get_threshold(wvfs, wvf_max .* 0.2; mintot=dsp_config.kwargs_pars.tx_mintot)
    t50 = get_threshold(wvfs, wvf_max .* 0.5; mintot=dsp_config.kwargs_pars.tx_mintot)
    t80 = get_threshold(wvfs, wvf_max .* 0.8; mintot=dsp_config.kwargs_pars.tx_mintot)
    t90 = get_threshold(wvfs, wvf_max .* 0.9; mintot=dsp_config.kwargs_pars.tx_mintot)
    t99 = get_threshold(wvfs, wvf_max .* 0.99; mintot=dsp_config.kwargs_pars.tx_mintot)
        
    drift_time = uconvert.(u"ns", t90 - t0)
    
    # get Q-drift parameter
    qdrift = get_qdrift(wvfs, t0, dsp_config.qdrift_int_length; pol_power=dsp_config.kwargs_pars.int_interpolation_order, sign_est_length=dsp_config.kwargs_pars.int_interpolation_length)
    
    # get LQ parameter
    lq  = get_qdrift(wvfs, t80, dsp_config.lq_int_length; pol_power=dsp_config.kwargs_pars.int_interpolation_order, sign_est_length=dsp_config.kwargs_pars.int_interpolation_length)
    
    # robust energy reconstruction with long, middle and short rise and flat-top times
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
    e_10410  = maximum.((uflt_10410.(wvfs)).signal)

    uflt_535 = TrapezoidalChargeFilter(5u"µs", 3u"µs")
    e_535  = maximum.((uflt_535.(wvfs)).signal)
    
    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")
    e_313  = maximum.((uflt_313.(wvfs)).signal)
    
    # signal estimator for precise energy reconstruction
    signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))
    
    # get trap energy of optimized rise and flat-top time
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
    e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))


    # get cusp energy of optimized rise and flat-top time
    uflt_cusp_rtft = CUSPChargeFilter(cusp_rt, cusp_ft, τ_cusp, flt_length_cusp, cusp_scale) 
    e_cusp = signal_estimator.(uflt_cusp_rtft.(wvfs), t50 .+ (flt_length_cusp /2))

    # get zac energy of optimized rise and flat-top time
    uflt_zac_rtft = ZACChargeFilter(zac_rt, zac_ft, τ_zac, flt_length_zac, zac_scale) 
    e_zac = signal_estimator.(uflt_zac_rtft.(wvfs), t50 .+ (flt_length_zac /2))


    # extract current with optimal SG filter length with second order polynominal and first derivative
    wvfs_sgflt_deriv = SavitzkyGolayFilter(sg_wl, dsp_config.sg_flt_degree, 1).(wvfs)
    a_sg = get_wvf_maximum.(wvfs_sgflt_deriv, leftendpoint(dsp_config.current_window), rightendpoint(dsp_config.current_window))
    
    a_60 = get_wvf_maximum.(SavitzkyGolayFilter(60u"ns", dsp_config.sg_flt_degree, 1).(wvfs), leftendpoint(dsp_config.current_window), rightendpoint(dsp_config.current_window))
    a_100 = get_wvf_maximum.(SavitzkyGolayFilter(100u"ns", dsp_config.sg_flt_degree, 1).(wvfs), leftendpoint(dsp_config.current_window), rightendpoint(dsp_config.current_window))
    a_raw = get_wvf_maximum.(DifferentiatorFilter(1).(wvfs), leftendpoint(dsp_config.current_window), rightendpoint(dsp_config.current_window))
    
    # get in-trace pile-up
    inTrace_pileUp = get_intracePileUp(wvfs_sgflt_deriv, dsp_config.inTraceCut_std_threshold, dsp_config.bl_window; mintot=dsp_config.kwargs_pars.intrace_mintot)
        
    # get position of current rise
    thres = maximum.(wvfs_sgflt_deriv.signal) .* 0.5
    # replace!(thres, zero(thres[1]) => one(thres[1]))
    
    t50_current = get_threshold(wvfs_sgflt_deriv, thres; mintot=dsp_config.kwargs_pars.tx_mintot)
    
    # output Table 
    # return 
    Table(blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        tailmean = pz_stats.mean, tailsigma = pz_stats.sigma, tailslope = pz_stats.slope, tailoffset = pz_stats.offset,
        t0 = t0, t10 = t10, t20 = t20, t50 = t50, t80 = t80, t90 = t90, t99 = t99,
        t50_current = t50_current, 
        drift_time = drift_time,
        tail_τ = tail_stats.τ, tail_mean = tail_stats.mean, tail_sigma = tail_stats.sigma,
        e_max = wvf_max, e_min = wvf_min,
        e_10410 = e_10410, e_535 = e_535, e_313 = e_313,
        e_trap = e_trap, e_cusp = e_cusp, e_zac = e_zac, 
        qdrift = qdrift, lq = lq,
        a_sg = a_sg, a_60 = a_60, a_100 = a_100, a_raw = a_raw,
        )
end 
export simple_dsp

"""
minimal version of `simple_dsp`  used to calculate quality cuts with peakfiles 
"""
function simple_dsp_qc(data::Q, dsp_config::DSPConfig; τ_pz::Quantity{T} = 0.0u"µs", pars_filter::PropDict = PropDict()) where {Q <: Table, T<:Real}
    wvfs = data.waveform

    #### get default from config 
    trap_rt, trap_ft = get_fltpars(pars_filter, :trap, dsp_config)

    ################## ACTUAL WAVEFORM FILTERING AND RECONSTRUCTION, ANALYSIS BEGINS HERE ##################
    bl_stats = signalstats.(wvfs, leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
   
    # tail analysis 
    tail_stats = tailstats.(wvfs, leftendpoint(dsp_config.tail_window), rightendpoint(dsp_config.tail_window))

    # deconvolute waveform: pole-zero correction. Use pre-defined tau from decay time analysis OR median decay time for all waveforms
    if τ_pz == 0.0u"µs"
        τ_pz = median(filter!(isfinite,tail_stats.τ))
    end
    deconv_flt = InvCRFilter(τ_pz)
    wvfs = deconv_flt.(wvfs)
  
    wvf_max = maximum.(wvfs.signal)

    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs, leftendpoint(dsp_config.tail_window), rightendpoint(dsp_config.tail_window))

    # t0  and t90 determination
    t0 = get_t0(wvfs, dsp_config.t0_threshold; flt_pars=dsp_config.kwargs_pars.t0_flt_pars, mintot=dsp_config.kwargs_pars.t0_mintot)
    t50 = get_threshold(wvfs, wvf_max .* 0.5; mintot=dsp_config.kwargs_pars.tx_mintot) 
    t90 = get_threshold(wvfs, wvf_max .* 0.9; mintot=dsp_config.kwargs_pars.tx_mintot)  
    drift_time = uconvert.(u"ns", t90 - t0)
    
    # energy from trap filter 
    signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
    e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))

    # output Table 
    Table(blmean = bl_stats.mean, 
        blslope = bl_stats.slope,
        tailmean = pz_stats.mean, 
        tailslope = pz_stats.slope,
        t0 = t0, 
        t90 = t90,
        tail_τ = tail_stats.τ,
        e_max = wvf_max, 
        e_trap = e_trap,
        drift_time = drift_time
        )
end 
export simple_dsp_qc

"""
    simple_dsp_pulser(data::Q, dsp_config::DSPConfig; τ_pz::Quantity{T} = 0.0u"µs", pars_filter::PropDict) where {Q <: Table, T<:Real}
minimal version of `simple_dsp` for pulser channel 
"""
function simple_dsp_pulser(data::Q, dsp_config::DSPConfig; τ_pz::Quantity{T} = 0.0u"µs", pars_filter::PropDict) where {Q <: Table, T<:Real}
    if !hasproperty(data, :pulser)
        error("data does not contain pulser waveforms")
    end 

    wvfs = data.pulser

    #### get filter parameters from config. if otimized values are not found; use default.  
    trap_rt, trap_ft = mvalue.(get_fltpars(pars_filter, :trap, dsp_config))

    ################## ACTUAL WAVEFORM FILTERING AND RECONSTRUCTION, ANALYSIS BEGINS HERE ##################
    bl_stats = signalstats.(wvfs, leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
   
    # tail analysis 
    tail_stats = tailstats.(wvfs, leftendpoint(dsp_config.tail_window), rightendpoint(dsp_config.tail_window))

    # get raw wvf maximum/minimum
    wvf_max = maximum.(wvfs.signal)
    wvf_min = minimum.(wvfs.signal)
    
    # deconvolute waveform: pole-zero correction. Use pre-defined tau from decay time analysis OR median decay time for all waveforms
    if τ_pz == 0.0u"µs"
        τ_pz = median(filter!(isfinite, tail_stats.τ))
    end
    deconv_flt = InvCRFilter(τ_pz)
    wvfs = deconv_flt.(wvfs)
  
    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs, leftendpoint(dsp_config.tail_window), rightendpoint(dsp_config.tail_window))

    # t0 determination
    t0 = get_t0(wvfs, dsp_config.t0_threshold; flt_pars=dsp_config.kwargs_pars.t0_flt_pars, mintot=dsp_config.kwargs_pars.t0_mintot)
    
    # if all waveforms are saturated set threshold to 1.0 to avoid numerical problems
    # replace!(wvf_max, zero(wvf_max[1]) => one(wvf_max[1]))
    
    # get threshold points in rise
    t10 = get_threshold(wvfs, wvf_max .* 0.1; mintot=dsp_config.kwargs_pars.tx_mintot)
    t20 = get_threshold(wvfs, wvf_max .* 0.2; mintot=dsp_config.kwargs_pars.tx_mintot)
    t50 = get_threshold(wvfs, wvf_max .* 0.5; mintot=dsp_config.kwargs_pars.tx_mintot)
    t80 = get_threshold(wvfs, wvf_max .* 0.8; mintot=dsp_config.kwargs_pars.tx_mintot)
    t90 = get_threshold(wvfs, wvf_max .* 0.9; mintot=dsp_config.kwargs_pars.tx_mintot)
    t99 = get_threshold(wvfs, wvf_max .* 0.99; mintot=dsp_config.kwargs_pars.tx_mintot)
        
    drift_time = uconvert.(u"ns", t90 - t0)
    
    # robust energy reconstruction with long, middle and short rise and flat-top times
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
    e_10410  = maximum.((uflt_10410.(wvfs)).signal)

    uflt_535 = TrapezoidalChargeFilter(5u"µs", 3u"µs")
    e_535  = maximum.((uflt_535.(wvfs)).signal)
    
    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")
    e_313  = maximum.((uflt_313.(wvfs)).signal)
    
    # signal estimator for precise energy reconstruction
    signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))
    
    # get trap energy of optimized rise and flat-top time
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
    e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))

    # output Table 
    # return 
    Table(pulser_blmean = bl_stats.mean, pulser_blsigma = bl_stats.sigma, pulser_blslope = bl_stats.slope, pulser_bloffset = bl_stats.offset, 
        pulser_tailmean = pz_stats.mean, pulser_tailsigma = pz_stats.sigma, pulser_tailslope = pz_stats.slope, pulser_tailoffset = pz_stats.offset,
        pulser_t0 = t0, pulser_t10 = t10, pulser_t20 = t20, pulser_t50 = t50, pulser_t80 = t80, pulser_t90 = t90, pulser_t99 = t99,
        pulser_drift_time = drift_time,
        pulser_tail_τ = tail_stats.τ, pulser_tail_mean = tail_stats.mean, pulser_tail_sigma = tail_stats.sigma,
        pulser_e_max = wvf_max, pulser_e_min = wvf_min,
        pulser_e_10410 = e_10410, pulser_e_535 = e_535, pulser_e_313 = e_313,
        pulser_e_trap = e_trap
        )
end 
export simple_dsp_pulser
