
"""
    charge of electron in Coulomb
"""
const electron_charge = 1.60217662e-19 # C

"""
GeBandgap(T::Real)
bandgap energy (eV) of germanium as a function of temperature (K)
"""
function GeBandgap(T::Real)
    return 0.744 - (4.774*1e-4 * T^2) / (T + 235)
end

"""
Ge_Energy_per_eholePair(T::Real)
energy (eV) required to create an electron-hole pair in germanium as a function of temperature (K)
"""
function GeEnergyer_eholePair(T::Real)
    return 2.2 * GeBandgap(T) + 1.99 * GeBandgap(T)^(3/2) * exp(4.75 * GeBandgap(T) / T)
end

"""
Ge_NumberChargeCarrier(E::Real, T::Real)
number of charge carriers created by an energy E (eV) in germanium at temperature T (K)
"""
function Ge_NumberChargeCarrier(E::Real, T::Real)
return E / Ge_Energy_per_eholePair(T)
end 

"""
    DAQ_ADC_to_V(ADC::Real, dynamicrange_V::Real, bits::Int)
Convert ADC-code from digitizer into voltage. The ADC-code is assumed to be in the range.
""" 
function DAQ_ADC_to_V(ADC::Real, dynamicrange_V::Real, bits::Int)
    return (dynamicrange_V / 2^bits) * ADC
end

"""
    pulser_ADC_to_electrons(ADC::Real, C_pulser::Real; bits::Int = 14, dynamicrange_V::Real = 2.0, gain::Real = 1.0)
Convert ADC-code from pulser into injected charge in pulser_ADC_to_electrons
"""
function pulser_ADC_to_electrons(ADC::Real, C_pulser::Real; bits::Int = 14,  dynamicrange_V::Real = 2.0, gain::Real = 1.0)
    V = DAQ_ADC_to_V(ADC, dynamicrange_V, bits)
    return  (V/gain * C_pulser) / electron_charge  # charge in electrons
end

"""
    pulser_ADC_to_keV(ADC::Real, C_pulser::Real; bits::Int = 14, dynamicrange_V::Real = 2.0, gain::Real = 1.0)
Convert ADC-code from pulser into injected energy in keV.
"""
function pulser_ADC_to_keV(ADC::Real, C_pulser::Real; bits::Int = 14,  dynamicrange_V::Real = 2.0, gain::Real = 1.0)
    return  pulser_ADC_to_electrons(ADC, C_pulser; bits = bits, dynamicrange_V = dynamicrange_V, gain = gain) * Ge_Energy_per_eholePair(90) / 1e3
end

