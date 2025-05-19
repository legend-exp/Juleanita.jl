# 💡 **Remarks**  
> This document provides a step-by-step guide for the initial setup of **Juleanita** on **NERSC**, using data from the **LBNL teststand** as an example.  
> If you plan to use a different set of teststand data, the overall logic remains the same, but you will need to adjust the file paths accordingly.

Working on [NERSC](https://www.nersc.gov/) is the most straighforward way to analyse LBNL teststand data, because:
1. Teststand data is stored on NERSC 
2. Metadata repository is already cloned into the correct location 
3. Configuration file for LegendDataManagement is already in place

## 1. NERSC account 
Get a NERSC account. For instructions see [Legend confluence](https://legend-exp.atlassian.net/wiki/spaces/LEGEND/pages/261750878/Computing+at+NERSC)

## 2. SSH connection
Once you have a NERSC account, you can log in via ssh: 
1. Create [MFA](https://docs.nersc.gov/connect/mfa/) (multi-factor authetification) to get OPT (one-time-password):
    * Once you have the account and MFA set up, you can login to NERSC:
    ``` bash
	   ssh <username>@perlmutter.nersc.gov
    ```
    replace `<username>` with your NERSC username.

2. Optional (but highly recommended): Download [sshproxy.sh](https://docs.nersc.gov/connect/mfa/#mfa-for-ssh-keys-sshproxy)
    * This creates ssh key that is valid 24 hours. No need to retype password+OTP everytime
	* usage:
     ``` bash
     ./sshproxy.sh -u <username>   
	```
3. Save connection in `~/.ssh/config` config to make connecting easier. Then login via `ssh perlmutter`
``` bash
Host perlmutter
    HostName perlmutter.nersc.gov
    User <username>
    IdentityFile ~/.ssh/nersc
    IdentitiesOnly yes
    ForwardAgent yes
```

## 3. VSCode 
Though not mandatory, I highly recommend to download and install the IDE [VSCode](https://code.visualstudio.com/). 
In addition, I'd recomment to install the following VSCode extentions: 
- `Julia`: language support
- `Remote - SSH`: allows you to connect to a remote, such as NERSC, inside VSCode.
- `GitHub Copilot` (requires GitHub Pro account): AI chat & coding help inside VSCode

To start a session on NERSC, open VSCode on your computer as usual. Then connect to NERSC via `Connect to Host` (bottom left) and select `perlmutter`. This will connect you to a NERSC login node. 

## 4. Julia on NERSC
### 4.1 Installation
You don't have to install Julia on NERSC yourself, because there are already several versions available, which are updated regularly.
See **/global/common/software/nersc/n9/julia/** for available versions. At the time of this writing the newest version is 1.10.4. 

### 4.2 Julia Executable 
When you start working with Julia on NERSC, some julia settings that have to be modified. You can access the julia settings in VSCode with 
**Settings -> Remote [SSH: permutter]-> Extensions -> julia**.

1. Julia Executable Path: this tells VSCode which julia installation you want to use (modify path if you want a different version). 
``` bash 
-> Julia: Executable Path =  /global/common/software/nersc/n9/julia/1.10.4/bin/julia
```
### 4.3 Start Julia session
To start a julia session, you can just type `julia` into the VSCode terminal window. Another way to start julia is via the Command Palette (Command + Shift + P) -> Julia: Start REPL. 

## 5. Julia package manager and Legend registry
### 5.1 Package manager
Julia has a built-in package manager. To enter the package manager inside  julia session type `]`
```julia
julia> ]
```
To get a list of packages in your julia environment type (from package manager):
```julia
pkg> st
```
To add a package to  your julia environment type (from package manager):
```julia
pkg> add <package-name>
```
### 5.2 Legend registry
To be able to add and use Legend-specific Julia packages, you need to add the Legend registry:  
```julia 
julia> include(download("https://raw.githubusercontent.com/legend-exp/legend-julia-tutorial/main/legend_julia_setup.jl"))
```
This needs to be done only once. For reference see [Julia Software Stack (Confluence)](https://legend-exp.atlassian.net/wiki/spaces/LEGEND/pages/494632973/Julia+Software+Stack#The-LEGEND-Julia-package-registry).

## 6. LBNL Teststand data 
### 6.1 Where is the data on NERSC? 
The LBNL teststand data is located on NERSC under this path: 
```bash
/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/
```

This folder contains
- **`config.json`**
    - this is a configuration file that defines where the raw data, intermediate and final results are saved
- **teststand-metadata**
    - Contains documentation and metadata for germanium detector and electronics
    - Contains configuration for processing (such as baseline window or gamma peak energies)
    - for more information see [git repo README](https://github.com/legend-exp/teststand-metadata)
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
            
### 6.2 How can I access the data?
To load and save data using the data management package `LegendDataManagement.jl` (a dependency of Juleanita), you need to define the environmental variable **`"LEGEND_DATA_CONFIG"`**. This variable should point to a `config.json` file in your data production. The config for the main LBNL production on NERSC is located here
```bash
/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json
```
This file tells LegendDataManagement where it can find the raw waveforms, dsp results, energy calibration function etc. 

You should define `"LEGEND_DATA_CONFIG"` in two places: 
1. your `.bashrc` as:

``` bash
export LEGEND_DATA_CONFIG="/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"
```
2. In your VSCode julia extension, you should also add the same config path to your VSCode [remote] settings. 
``` json
"terminal.integrated.env.linux": {
            "LEGEND_DATA_CONFIG": "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"
        },
```

##  Now you are ready to start using Juleanita with your teststand data. 