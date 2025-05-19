# read_data
A basic introduction on how to read data using LegendDataManagement

First you have to load the required packages

```` reading_data
using LegendDataManagement
using LegendHDF5IO
using TypedTables
````

Now you define which data set you want to read.
For more information on the datasets, check out the documentation in [teststand-metadata](https://github.com/legend-exp/teststand-metadata/tree/main/ppc01/doc)

```` reading_data
category = DataCategory(:cal)
period = DataPeriod(3)
run = DataRun(1)
channel = ChannelId(1) # channel 1 => PPC01 in LBNL-70-141
````

Create LegendData object

```` reading_data
asic = LegendData(:ppc01)
````

## Read the "raw" tier (waveforms)

```` reading_data
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
data_raw = Table(read_ldata(asic, DataTier(:raw), filekeys, channel))
waveforms = data_raw.waveform
````

## Read the "dsp" tier

````  reading_data
data_dsp = Table(read_ldata(asic, DataTier(:jldsp), filekeys, channel))
````

this data contains many variables, like the different energy estimates such as

```` reading_data
data_dsp.e_trap
````

## Read the "hit" tier

````reading_data
data_hit = Table(read_ldata(asic, DataTier(:jlhit), filekeys, channel))
````

this data contains the hit information, such as the calibrated energyes and quality cuts flags

## Read the "pars"
Saved parameters, can be access with the following synthax:

````reading_data
pars = asic.par[category].rpars
````

here are some example to lad specific parameters:

```` reading_data
pars_pole_zero = pars.pz[period, run, channel]
pars_energy = pars.ecal[period, run, channel].e_trap
gamma_peaks = [pars_energy.fit[peak].fwhm for peak in [:Co60a, :Co60b]]
qc_flag =  pars.qc[period, run, channel].wvf_keep.all
````

