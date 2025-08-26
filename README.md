# PeriRxn

Python scripts for creating Gaussian input files for analyzing **NICSzz** and **normalized multicentered bond orders (nMCBO)** to characterize pericyclic reactions.  

The package also provides utilities for figure generation, such as the examples shown below:

![mNICS + mNMCBO Plot](https://github.com/WillDeSnoo/PeriRxn/blob/02_01_2025/exampleCope/mNICS_mNMCBO_IRC_Plot.png)  
![NICS IRC Scan Plot](https://github.com/WillDeSnoo/PeriRxn/blob/02_01_2025/exampleCope/NICS_IRC_Scan_Plot.png)

---

## Features
- Generate Gaussian input files for NICSzz and nMCBO calculations.  
- Analyze reaction pathways and extract bonding indicators.  
- Plot NICS and nMCBO values along intrinsic reaction coordinates (IRCs).  
- Visualize aromatic/antiaromatic character during pericyclic reactions.  

---

## Installation

Clone the repository:
```bash
git clone https://github.com/WillDeSnoo/PeriRxn.git
cd PeriRxn
```
Make the main scripts executable:
```
chmod u+x perirxn_2.0.py
chmod u+x peridata_ShieldTensor.py
```
Add the scripts to your $PATH (e.g. in ~/.bashrc or ~/.zshrc):
```
export PATH="$PATH:/path/to/PeriRxn"
```
---
## Requirements

Python 3.8+

Multiwfn
 (required for the nMCBO portion of the analysis)

## Usage
### Step 1: Generate Gaussian Input Files

From a Gaussian IRC log file, run:
```
perirxn_2.0.py irc.log
```

Notes:

The IRC must be bidirectional.

If the level of theory is not detected, you will be prompted to enter it.

You will also be asked to specify indices of aromatic atoms:
```
Please type the indices of the aromatic carbons (i.e. 1,2,3,4...):
If Diels–Alder type reaction, enter the indices of aromatic atoms of two fragments separated by a space (i.e. 1,2,3,4 5,6):
```

This creates a directory called peri_irc/ containing input files with regularly spaced ghost atoms for NICSzz sampling. Run these input files with Gaussian.


### Step 2: Analyze and Plot Results

After running the Gaussian input files, move into the peri_irc/ directory containing the .log files and run:
```
peridata_ShieldTensor.py irc.log
```
This script will:

Extract NICS and nMCBO data

Generate plots (as shown above)

Save results to nics_mcbo_data.xlsx for further analysis
