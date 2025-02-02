#!/usr/bin/env python3

#Imports
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import os

def read_xlsx_sheet(fname, id):
    xl = pd.ExcelFile(fname)
    sheet_names=xl.sheet_names
    lines = pd.read_excel(fname, sheet_name = id)
    return lines, sheet_names

def get_data_from_excel(fname, sheet_num=0,row0=0):
    sheet_names = []
    df, sheet_names = read_xlsx_sheet(fname, sheet_num)
    print("SHEET NAMES:", sheet_names)
    df_array = df.to_numpy()
    
    irc_step_array = df_array[row0:,0]
    irc_val_array = df_array[row0:,1]
    energy_array=df_array[row0:,3]
    nicszz = df_array[row0:,4]
    mcbo = df_array[row0:, 5]
    return irc_val_array, energy_array, nicszz, mcbo, irc_step_array

def get_heatmap_data_from_excel(fname, sheet_num=1):
    df, sheet_names = read_xlsx_sheet(fname, sheet_num)

def plot_data1(irc_val_array, energy_array, nicszz, mcbo=None, irc_step_array=None, savedir='/periplots',outname="1",verbose=False):
    
    #irc_val_array = irc_step_array #Uncomment this line to use IRC step as x-axis instead of IRC val
    
    nicszz_mask = []
    mcbo_mask = []
    
    for i,val in enumerate(nicszz):
        if not val:
            nicszz[i]=np.nan
            nicszz_mask.append(False)
        else:
            nicszz_mask.append(True)
    
    for i,val in enumerate(mcbo):
        if not val:
            mcbo[i]=np.nan
            mcbo_mask.append(False)
 
        else:
            mcbo_mask.append(True)
    
    fig, ax = plt.subplots(figsize=[6,4])
    fig.subplots_adjust(right=0.75)
    
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    
    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines['right'].set_position(("axes", 1.2))
    
    p1, = ax.plot(irc_val_array, energy_array, "-", color = 'blue', label="Energy")
    if verbose:
        print("TEST")
        print("nicszz_mask:", nicszz_mask)
        print(irc_val_array)
        print(nicszz[nicszz_mask])
    p2, = twin1.plot(irc_val_array[nicszz_mask], nicszz[nicszz_mask], "-", color = 'red', label=r'$\int$NICS$_{zz}$ (ppm)')
    if mcbo.any:
        p3, = twin2.plot(irc_val_array[mcbo_mask], mcbo[mcbo_mask], "-", color = 'green', label="nMCBO")
    
    #ax.set_title('{}'.format(sheet_names[sheet_num]))
    ax.set_xlim(min(irc_val_array),max(irc_val_array))
    ax.set_ylim(min(energy_array)-0.05*max(abs(energy_array)), max(energy_array)+0.05*max(abs(energy_array)))
    twin1.set_ylim(min(nicszz),max(nicszz))
    
    if mcbo.any:
        twin2.set_ylim(min(mcbo),max(mcbo))
    
    #twin1.set_ylim(-20, 20) # ADJUST THESE AS NEEDED
    #twin2.set_ylim(-0.8, 0.8) # ADJUST THESE AS NEEDED
    
    twin1.set_ylim(min(nicszz)-0.05*max(abs(nicszz)), max(nicszz)+0.05*max(abs(nicszz))) # ADJUST THESE AS NEEDED
    twin2.set_ylim(min(mcbo)-0.05*min(abs(mcbo)), max(mcbo)+0.05*max(abs(mcbo))) # ADJUST THESE AS NEEDED
    
    
    ax.set_xlabel("Intrinsic Reaction Coordinate")
    ax.set_ylabel("Energy (kcal/mol)")
    twin1.set_ylabel(r'$\overline{NICS}$$_{zz}$ (ppm)')
    if mcbo.any:
        twin2.set_ylabel("nMCBO")
    
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    if mcbo.any:
        twin2.yaxis.label.set_color(p3.get_color())
    
    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='x', colors = 'black', **tkw)
    ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    if mcbo.any:
        twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    plt.show()
    #plt.savefig(savedir+'/Periplot_{}'.format(outname), dpi=300)

def plot_data(irc_val_array, energy_array, nicszz, mcbo=None, irc_step_array=None, savedir='/periplots',outname="1",verbose=False):
    print('Plotting data...')
    #irc_val_array = irc_step_array #Uncomment this line to use IRC step as x-axis instead of IRC val
    
    nicszz_mask = []
    mcbo_mask = []

    # Create masked data, to plot only data with values
    for i,val in enumerate(nicszz):
        if not val or val != val:
            nicszz[i]=np.nan
            nicszz_mask.append(False)
        else:
            nicszz_mask.append(True)
    
    for i,val in enumerate(mcbo):
        if not val or val != val:
            mcbo[i]=np.nan
            mcbo_mask.append(False)
 
        else:
            mcbo_mask.append(True)
    
    nicszz_mask=mcbo_mask
    # Create plot and axes with dimensions.
    fig, ax = plt.subplots(figsize=[6,4])
    fig.subplots_adjust(right=0.75)
    
    # Create twin axes for nics and mcbo
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    
    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines['right'].set_position(("axes", 1.2))
    
    p1, = ax.plot(irc_val_array, energy_array, "-", color = 'blue', label="Energy")
    if verbose:
        print("verbose")
        print("TEST")
        print("nicszz_mask:", nicszz_mask)
        print("nmcbo_mask", mcbo_mask)
        print(irc_val_array)
        print("NICS values:")
    
    p2, = twin1.plot(irc_val_array[nicszz_mask], nicszz[nicszz_mask], "-", color = 'red', label=r'$\overline{\text{NICS}_{zz}} (ppm')
    #p2, = twin1.plot(irc_val_array[nicszz_mask], nicszz[nicszz_mask], "-", color = 'red', label=r'$\int$NICS_{zz}$$ (ppm)')
    if mcbo.any:
        p3, = twin2.plot(irc_val_array[mcbo_mask], mcbo[mcbo_mask], "-", color = 'green', label="nMCBO")
    
    #ax.set_title('{}'.format(sheet_names[sheet_num]))
    ax.set_xlim(min(irc_val_array),max(irc_val_array))
    max_abs_energy_array=max(abs(np.asarray(energy_array)))
    ax.set_ylim(min(energy_array)-0.05*max_abs_energy_array, max(energy_array)+0.05*max_abs_energy_array)
    
    twin1.set_ylim(min(nicszz),max(nicszz))
    if mcbo.any:
        twin2.set_ylim(min(mcbo),max(mcbo))
    
    twin1.set_ylim(-30, 30) # ADJUST THESE AS NEEDED
    twin2.set_ylim(-1.0, 1.0) # ADJUST THESE AS NEEDED
    max_abs_nicszz=max(abs(nicszz))
    #twin1.set_ylim(-max_abs_nicszz-0.05*max_abs_nicszz, max_abs_nicszz+0.05*max_abs_nicszz) # ADJUST THESE AS NEEDED
    #twin2.set_ylim(min(mcbo)-0.05*min(abs(mcbo)), max(mcbo)+0.05*max(abs(mcbo))) # ADJUST THESE AS NEEDED
    
    
    ax.set_xlabel("Intrinsic Reaction Coordinate")
    ax.set_ylabel("Energy (kcal/mol)")
    twin1.set_ylabel(r'$\overline{NICS_{zz}}$ (ppm)')
    if mcbo.any:
        twin2.set_ylabel("nMCBO")
    
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    if mcbo.any:
        twin2.yaxis.label.set_color(p3.get_color())
    
    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='x', colors = 'black', **tkw)
    ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    if mcbo.any:
        twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    plt.savefig("./mNICS_mNMCBO_IRC_Plot")
    plt.show()

    #plt.savefig(savedir+'/Periplot_{}'.format(outname), dpi=300)

def plotNICSHM(irc_val_array, energy_array, allnicszz, irc_step, bqarray=np.arange(-3,3.1,0.1)):
    #x=irc_val_array
    x=irc_step
    y=bqarray
    Z=np.transpose(allnicszz)
    X, Y = np.meshgrid(x, y)
    # Create a figure and plot the paraboloid
    plt.figure()
    plt.contourf(X, Y, Z, levels=20, cmap='viridis')
    plt.colorbar(label='NICS$_{zz}$ (ppm)')
    plt.xlabel('Intrinsic Reaction Coordinate')
    plt.ylabel('Z distance (Å)')
    #plt.title('NICSzz IRC Scan Plot')
    plt.savefig("./NICS_IRC_Scan_Plot")
    plt.show()
    
def plotFromExcel(filename):
    excel_data = get_data_from_excel(filename)
    #plotNICSHM(np.asarray(column(energies, 0)), energiesKcalMol,np.asarray(full_bq_array),irc_step)

    #irc_val_array, energy_array, nicszz, mcbo, irc_step_array = get_data_from_excel(filename)
    plot_data(*excel_data)
    
def main(argv):
    plotFromExcel(argv)
    
if __name__ == "__main__":
    main(sys.argv[1])
    
    
