# MakeSumPeakHistograms
Creates energy-angle plots for nuclear two-photon decay. Time random background
subtraction occurs while the histograms are being filled and a second event
mixed (depth 10) matrix is created for normalization procedures.

## Table of Contents
  * [Installation](#installation)
  * [Running MakeSumPeakHistograms](#running-makesumpeakhistograms)
    + [Setting Energy Gates](#settingenergygates)
  * [Helper scripts](#helper-scripts)
    + [MakeSumPeakHistograms.sh](#makesumpeakhistogramssh)

# Installation
0. Requires GRSISort 4.X.X
1. Get the code, either via tarball or from github
```
git clone https://github.com/cnatzke/SumPeakAnalysis.git
```
2. Build program using standard cmake process, e.g.
```
mkdir myBuild && cd myBuild
cmake ..
make
make install # if you would like to install the programs elsewhere.
```
3. Do science.

# Running MakeSumPeakHistograms
The general form of input is:
```
./SumPeakHistograms analysis_tree [analysis_tree_2 ... ] calibration_file linear_parameter_file
```

##### Parameters
```
analysis_tree           ROOT file(s) containing analysis tree to process (must end with .root)
calibration_file        GRIFFIN calibration file (must end with .cal)
linear_parameter_file   File containing secondary linear calibration coefficients (must end with .txt)
```

##### Outputs
```
secondary_calibrated_histograms.root   ROOT file containing energy-angle plots
```

## Setting Energy Gates
The energies of interest are defined in ```HistogramManager.cpp``` file vector ```energy_gates_vec```. To add an energy plot simply add the energy to the vector and recompile the program.


# Helper scripts
Included is a helper script that makes building histograms easier.

### MakeSumPeakHistograms.sh
This script takes a run and a batch size.
```
./MakeSumPeakHistograms.sh run_number batch_size
```
where
```
run_number  Run number of interest
batch_size  Number of subruns to be chained together
```
The script chains subruns together to process, but depending on the size of the run  this exceeds the available system memory. The batch size controls how many subruns will be chained together in each step. For example if you have 17 subruns (0-16) and a batch size of 3, the script will pass a total of 6 batches to the histogram programs like so
```
batch   subruns
---------------
1       000-002
2       003-005
3       006-008
4       009-011
5       012-014
6       015-016
```
