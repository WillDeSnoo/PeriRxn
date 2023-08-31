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

def get_data_from_excel(fname, sheet_num=0):
    sheet_names = []
    df, sheet_names = read_xlsx_sheet(fname, sheet_num)
    print("SHEET NAMES:", sheet_names)
    df_array = df.to_numpy()
    
    headers = df_array[2,:]
    irc_step_array = df_array[3:,0]
    irc_val_array = df_array[3:,1]
    energy_array=df_array[3:,3]
    nicszz = df_array[3:,4]
    mcbo = df_array[3:, 5]
    allnicszz=df_array[3:,6:66]
    return irc_val_array, energy_array, nicszz, mcbo, irc_step_array

def plot_data(irc_val_array, energy_array, nicszz, mcbo=None, irc_step_array=None, savedir='/periplots',outname="1"):
    
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
    
    masked_nicszz = np.ma.masked_where(nicszz != nicszz, nicszz) # If entry is nan, it will be masked
    
    fig, ax = plt.subplots(figsize=[6,4])
    fig.subplots_adjust(right=0.75)
    
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    
    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines['right'].set_position(("axes", 1.2))
    
    p1, = ax.plot(irc_val_array, energy_array, "-", color = 'blue', label="Energy")
    print("TEST")
    print("nicszz_mask:", nicszz_mask)
    print(irc_val_array)
    print(nicszz[nicszz_mask])
    p2, = twin1.plot(irc_val_array[nicszz_mask], nicszz[nicszz_mask], "-", color = 'red', label=r'$\int$NICS$_{zz}$ (ppm)')
    if mcbo:
        p3, = twin2.plot(irc_val_array[mcbo_mask], mcbo[mcbo_mask], "-", color = 'green', label="nMCBO")
    
    #ax.set_title('{}'.format(sheet_names[sheet_num]))
    ax.set_xlim(min(irc_val_array),max(irc_val_array))
    
    
    
    ax.set_ylim(min(energy_array), max(energy_array))
    print(nicszz)
    twin1.set_ylim(min(nicszz),max(nicszz))
    
    if mcbo:
        twin2.set_ylim(min(mcbo),max(mcbo))
    
    #twin1.set_ylim(-20, 20) # ADJUST THESE AS NEEDED
    #twin2.set_ylim(-0.8, 0.8) # ADJUST THESE AS NEEDED
    
    twin1.set_ylim(min(nicszz)-0.05*max(abs(nicszz)), max(nicszz)+0.05*max(abs(nicszz))) # ADJUST THESE AS NEEDED
    twin2.set_ylim(-20, 20) # ADJUST THESE AS NEEDED
    
    
    ax.set_xlabel("Intrinsic Reaction Coordinate")
    ax.set_ylabel("Energy (kcal/mol)")
    twin1.set_ylabel(r'$\overline{NICS}$$_{zz}$ (ppm)')
    if mcbo:
        twin2.set_ylabel("nMCBO")
    
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    if mcbo:
        twin2.yaxis.label.set_color(p3.get_color())
    
    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='x', colors = 'black', **tkw)
    ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    if mcbo:
        twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    '''
    try:  # Tries creating new directory to place .com files in
        os.makedirs(savedir)  # Making directory for files
    except FileExistsError:
        pass
    '''
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
    
    