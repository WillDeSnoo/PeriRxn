# PeriRxn
Python scripts for creating gaussian input files for analyzing NICSzz and normalized multicentered bond orders (nMCBO) to characterize pericyclic reactions. Create figures like the one shown below:
![alt text](https://github.com/WillDeSnoo/PeriRxn/blob/02_01_2025/exampleCope/mNICS_mNMCBO_IRC_Plot.png)
![alt text](https://github.com/WillDeSnoo/PeriRxn/blob/02_01_2025/exampleCope/NICS_IRC_Scan_Plot.png)

---

## Features
- Generate Gaussian input files for NICSzz and nMCBO calculations.  
- Analyze reaction pathways and extract bonding indicators.  
- Plot NICS and nMCBO values along intrinsic reaction coordinates (IRCs).  
- Tools for visualizing aromatic/antiaromatic character during pericyclic reactions.  

---

## Installation

Clone the repository:
```bash
git clone https://github.com/WillDeSnoo/PeriRxn.git
cd PeriRxn

First make the following scripts executable
chmod u+x perirxn_2.0.py
chmod u+x peridata_ShieldTensor.py

Then add the scripts to path.
export PATH="$PATH:/path/to/PeriRxn/scripts"


---

## Usage

Create input files from a gaussian irc log file with perirxn_2.0.py. Note the irc should be bidirectional. Answer the respective prompts.

perirxn_2.0.py irc.log

1,2,3,4,5,6

This will create a directory called peri_irc with input files with reguarly spaced ghost atoms to sample NICSzz. Run these input files with Gaussian.

After running the newly generated inputfiles, we can analyze and plot the data with the peridata_ShieldTensor.py script. Enter the peri_irc directory with the log files, execute the script and follow the prompts.

peridata_ShieldTensor.py 

Generate plots as shown above, and save the data to nics_mcbo_data.xlsx file for further analysis.
