import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

print("Running periplot.py...")

class PeriData:
    def __init__(self, headers, irc_val_array, energy_array, nicszz, mcbo, irc_step_array, allnicszz):
        self.headers = np.array(headers)
        self.irc_val_array = np.array(irc_val_array)
        self.energy_array = np.array(energy_array)
        self.nicszz = np.array(nicszz)
        self.mcbo = np.array(mcbo)
        self.irc_step_array = np.array(irc_step_array)
        self.allnicszz = np.array(allnicszz)
    
    def plot_data(self):
        nicszz_mask=np.isfinite(self.nicszz)
        mcbo_mask=np.isfinite(self.mcbo)
    
        fig, ax = plt.subplots(figsize=[6,4])
        fig.subplots_adjust(right=0.75)
        
        twin1 = ax.twinx()
        
        # Offset the right spine of twin2.  The ticks and label have already been
        # placed on the right by twinx above.
        
        
        p1, = ax.plot(self.irc_val_array, self.energy_array, "-", color = 'blue', label="Energy")
        p2, = twin1.plot(self.irc_val_array[nicszz_mask], self.nicszz[nicszz_mask], "-", color = 'red', label=r'$\int$NICS$_{zz}$ (ppm)')
        tkw = dict(size=4, width=1.5)
        if not np.isnan(self.mcbo).all():
            print("nMCBO found")
            twin2 = ax.twinx()
            twin2.spines['right'].set_position(("axes", 1.2))
            p3, = twin2.plot(self.irc_val_array[mcbo_mask], self.mcbo[mcbo_mask], "-", color = 'green', label="nMCBO")
            twin2.set_ylim(-1.0,1.0)
            twin2.set_ylabel("nMCBO")
            twin2.yaxis.label.set_color(p3.get_color())
            twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
        else:
            print("No nMCBO detected")
        
        #ax.set_title('{}'.format(sheet_names[sheet_num]))
        ax.set_xlim(min(self.irc_val_array),max(self.irc_val_array))
        ax.set_ylim(min(self.energy_array), max(self.energy_array))
        
        twin1.set_ylim(-20, 20) # ADJUST THESE AS NEEDED
        
        ax.set_xlabel("Intrinsic Reaction Coordinate")
        ax.set_ylabel("Energy (kcal/mol)")
        twin1.set_ylabel(r'$\overline{NICS}$$_{zz}$ (ppm)')
        
        ax.yaxis.label.set_color(p1.get_color())
        twin1.yaxis.label.set_color(p2.get_color())
        
        ax.tick_params(axis='x', colors = 'black', **tkw)
        ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
        twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        '''
        try:  # Tries creating new directory to place .com files in
            os.makedirs(savedir)  # Making directory for files
        except FileExistsError:
            pass
        '''
        plt.show()
        #plt.savefig(savedir+'/Periplot_{}'.format(outname), dpi=300)
        
    def plotNICSHM(self):
        print(f'allNICSzz {self.allnicszz}')
        validIndicies=np.argwhere(np.isfinite(self.nicszz)).flatten()
        x=self.irc_val_array[validIndicies]
        print(f'valInd: {validIndicies}')
        y=np.arange(-3,3.1,0.1) #Change this if bq placement was altered from default
        print(f'irc_val_array: {x}')
        noNanNICS=self.allnicszz[~np.isnan(self.allnicszz)].reshape(len(x),len(y))
        print(f'noNanNICS:{noNanNICS}')
        Z=np.transpose(noNanNICS)
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
        
def read_xlsx_sheet(fname, id):
    xl = pd.ExcelFile(fname)
    sheet_names=xl.sheet_names
    lines = pd.read_excel(fname, sheet_name = id)
    return lines, sheet_names

def get_data_from_excel(fname, sheet_num=0):
    sheet_names = []
    df, sheet_names = read_xlsx_sheet(fname, sheet_num)
    df_array = df.to_numpy()
    
    headers = df_array[0,:]
    irc_step_array = df_array[0:,0]
    irc_val_array = df_array[0:,1]
    energy_array=df_array[0:,3]
    nicszz = df_array[0:,4]
    mcbo = df_array[0:, 5]
    allnicszz = df_array[0:,6:67]
    return headers, irc_val_array, energy_array, nicszz, mcbo, irc_step_array, allnicszz
    
def plot_data(irc_val_array, energy_array, nicszz, mcbo=None, irc_step_array=None, savedir='/periplots',outfilename="nMCBO_NICSzz_1D_plot"):
    #irc_val_array = irc_step_array #Uncomment this line to use IRC step as x-axis instead of IRC val
    nicszz_mask=np.isfinite(nicszz)
    mcbo_mask=np.isfinite(mcbo)

    fig, ax = plt.subplots(figsize=[6,4])
    fig.subplots_adjust(right=0.75)
    
    twin1 = ax.twinx()
    
    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    
    
    p1, = ax.plot(irc_val_array, energy_array, "-", color = 'blue', label="Energy")
    p2, = twin1.plot(irc_val_array[nicszz_mask], nicszz[nicszz_mask], "-", color = 'red', label=r'$\int$NICS$_{zz}$ (ppm)')
    tkw = dict(size=4, width=1.5)
    if not np.isnan(mcbo).all():
        print("nMCBO found")
        twin2 = ax.twinx()
        twin2.spines['right'].set_position(("axes", 1.2))
        p3, = twin2.plot(irc_val_array[mcbo_mask], mcbo[mcbo_mask], "-", color = 'green', label="nMCBO")
        twin2.set_ylim(-0.8,0.8)
        twin2.set_ylabel("nMCBO")
        twin2.yaxis.label.set_color(p3.get_color())
        twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    else:
        print("No nMCBO detected")
    
    #ax.set_title('{}'.format(sheet_names[sheet_num]))
    ax.set_xlim(min(irc_val_array),max(irc_val_array))
    ax.set_ylim(min(energy_array), max(energy_array))
    
    twin1.set_ylim(-20, 20) # ADJUST THESE AS NEEDED
    
    ax.set_xlabel("Intrinsic Reaction Coordinate")
    ax.set_ylabel("Energy (kcal/mol)")
    twin1.set_ylabel(r'$\overline{NICS}$$_{zz}$ (ppm)')
    
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    
    ax.tick_params(axis='x', colors = 'black', **tkw)
    ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    if outfilename:
        plt.savefig(outfilename, dpi=300)
    plt.show()
    
    
def plotNICSHM(irc_val_array, energy_array, allnicszz, irc_step, bqarray=np.arange(-3,3.1,0.1),outfilename='NICSzz_2D_Plot'):
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
    if outfilename:
        plt.savefig(outfilename)
    plt.show()
    
def main():
    parser = argparse.ArgumentParser(description="A script reads data from .xlsx and plots mean NICS or a 2d Heat map of NICS data along an IRC.")

    # Positional argument for the filename (required).
    parser.add_argument("filename", help="The name of the file containing data (.xlsx). Use peridata to generate this file")

    # Optional argument for an integer value."
    parser.add_argument("-plotType", "-pt", type=int,nargs='?', const=1, default=1, help="plotType: if 1: Plot mNICS and nMCBO. 2: Plot 2D Plot of NICS (x time, y is displacement above ring).")
    
    parser.add_argument("-save", "-s", type=str, nargs='?', const=None, default=None, help="Output file name to save plot")
    args = parser.parse_args()
    #End of argparsing
    
    #Get data from parsing
    plotDat=get_data_from_excel(args.filename)
    #Save data into a Class
    periData=PeriData(*plotDat)
    #Plot Data
    if args.plotType==1:
        periData.plot_data()
    elif args.plotType==2:
        periData.plotNICSHM()
    else:
        print("Invalid plotType option. Try periplot.py -h for help.")
 
if __name__ == "__main__":
    main()
    
    
    