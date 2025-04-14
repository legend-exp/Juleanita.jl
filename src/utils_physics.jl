
### Constants
"""
    electron_charge

charge of electron in Coulomb: ``e = 1.60217662 10^{-19}`` C
"""
const electron_charge = 1.60217662e-19
export electron_charge

### Germanium specifics 
"""
    GeBandgap(T::Real)
bandgap energy (eV) of germanium as a function of temperature (K)

```math 
E_\\textrm{gap}(T) = 0.744 - \\frac{4.774 \\cdot 10^{-4} \\cdot T^2}{T + 235}
```
"""
function GeBandgap(T::Real)
    return 0.744 - (4.774*1e-4 * T^2) / (T + 235)
end

"""
    Ge_Energy_per_eholePair(T::Real)
energy (eV) required to create an electron-hole pair in germanium as a function of temperature (K)

```math 
E(T)  = 2.2 \\cdot E_\\textrm{gap}(T) + 1.99 \\cdot E_\\textrm{gap}(T)^{3/2} \\cdot \\exp(4.75 \\cdot \\frac{E_\\textrm{gap}(T)}{T})
```
"""
function Ge_Energy_per_eholePair(T::Real)
    return 2.2 * GeBandgap(T) + 1.99 * GeBandgap(T)^(3/2) * exp(4.75 * GeBandgap(T) / T)
end

"""
    Ge_NumberChargeCarrier(E::Real, T::Real)
number of charge carriers created by an energy E (eV) in germanium at temperature T (K)
"""
function Ge_NumberChargeCarrier(E::Real, T::Real)
return E / Ge_Energy_per_eholePair(T)
end 

###  Conversion util functions (ADC - Volts - electrons - keV)
"""
    _ADC_to_V(ADC::Real, dynamicrange_V::Real, bits::Int)
Convert ADC-code from digitizer into voltage. The ADC-code is assumed to be in the range [0, 2^bits] corresponding to voltages within the dynamic range.
INPUTS:
- `ADC`: ADC-code from digitizer
- `dynamicrange_V`: dynamic range of the DAQ in Volts
- `bits`: number of bits of the ADC
""" 
function _ADC_to_V(ADC::Real, dynamicrange_V::Real, bits::Int)
    return (dynamicrange_V / 2^bits) * ADC
end

"""
    _ADC_to_electrons(ADC::Real, capacitance_F::Real; bits::Int = 14, dynamicrange_V::Real = 2.0, gain::Real = 1.0)
Convert ADC-code from pulser into injected charge in pulser_ADC_to_electrons
Inputs:
- `ADC``: ADC-code from digitizer  
- `capacitance_F`: capacitance of the system (could be pulser, detector, or combination) in Farad
- `bits`: number of bits of the ADC
- `dynamicrange_V`: dynamic range of the DAQ in Volts
- `gain`: gain of the system
"""
function _ADC_to_electrons(ADC::Real, capacitance_F::Real; bits::Int = 14,  dynamicrange_V::Real = 2.0, gain::Real = 1.0)
    V = _ADC_to_V(ADC, dynamicrange_V, bits)
    return  (V/gain * capacitance_F) / electron_charge  # charge in electrons
end

"""
    pulser_ADC_to_keV(ADC::Real, capacitance_F::Real; bits::Int = 14, dynamicrange_V::Real = 2.0, gain::Real = 1.0)
Convert ADC-code from pulser into injected energy in keV.
Inputs: 
- `ADC`: ADC-code from digitizer
- `capacitance_F`: capacitance of the system (could be pulser, detector, or combination) in Farad
- `bits`: number of bits of the DAQ
- `dynamicrange_V`: dynamic range of the ADC in Volts
- `gain`: gain of the system
"""
function _ADC_to_keV(ADC::Real, capacitance_F::Real; bits::Int = 14,  dynamicrange_V::Real = 2.0, gain::Real = 1.0)
    return _ADC_to_electrons(ADC, capacitance_F; bits = bits, dynamicrange_V = dynamicrange_V, gain = gain) * Ge_Energy_per_eholePair(90) / 1e3
end

"""
    V_to_electrons(Voltage_V::Real, capacitance_F::Real; gain::Real = 1.0)
Convert voltage into a charge in electrons based on the capacitance of the system (could be pulser, detector, or combination).
"""
function V_to_electrons(Voltage_V::Real, capacitance_F::Real; gain::Real = 1.0)
    return (Voltage_V/gain * capacitance_F) / electron_charge  # charge in electrons
end