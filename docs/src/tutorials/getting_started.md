# Getting started - Data

## Where is the data located (on NERSC)?
The LBNL teststand data is located on NERSC under this path: 
```
/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/
```

This folder contains
- **`config.json`**
    - this is a configuration file that defines where the raw data, intermediate and final results are saved
- **teststand-metadata**
    - Contains documentation and metadata for germanium detector and electronics
    - Contains configuration for processing (such as baseline window or gamma peak energies)
    - for more information see [git repo README](https://github.com/LisaSchlueter/teststand-metadata)
- **generated** 
    - this is where all the data lies
        - `jlpar`: results from processors, for example: decay time or energy calibration function
        - `jlplt`: processing plots
        - `jlreport`: space for processing reports. Not used at the moment. 
        - `tier/`: 
            - `jldsp/`: dsp results
            - `jlhit/`: hit files (calibrated spectra)
            - `jlpeaks/` peak files
            - `raw_csv/` waveforms in csv format (from DAQ)
            - `raw/` waveforms better compressed in `lh5` format (this is what LegendDataManagement needs)
            
### How can I access the data?
To load and save data using the data management package `LegendDataManagement.jl` (a dependecy of Juleanity), you need to define the environmental variable `"LEGEND_DATA_CONFIG"`. This variable should point to a `config.json` file in your data production. 

For example, on NERSC for the main data production you should define it in your `.bashrc` as:

``` bash
export LEGEND_DATA_CONFIG="/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"
```

If you use VSCode with the julia extension, you should also add the same config path to your VSCode [remote] settings. 
``` json
"terminal.integrated.env.linux": {
            "LEGEND_DATA_CONFIG": "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"
        },
```

## Legend packages
Juleanita already has many relevant Legend-specific packages as its dependency. If you want to load additionaly Legend packages, you have to add the Legend-Julia registry by executing the following command in a julia session:
``` julia
include(download("https://raw.githubusercontent.com/legend-exp/legend-julia-tutorial/main/legend_julia_setup.jl"))
```

This needs to be done only once. For reference see [confluence](https://legend-exp.atlassian.net/wiki/spaces/LEGEND/pages/494632973/Julia+Software+Stack#The-LEGEND-Julia-package-registry).
