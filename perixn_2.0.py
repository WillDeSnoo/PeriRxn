#!/usr/bin/env python3

#Imports
import numpy as np
import periodictable as pt
import sys
import os
import matplotlib.pyplot as plt
import re
import shutil

# CHANGE FOLLOWING PARAMETERS ACCORDINGLY
# -------------------------------------------------------------------------- #

# Default infilename, will be overidden by commandline input
infilename="/Users/willdesnoo/scripts/test/test_irc.log"

# Full default, set to False to prompt user to check/override default parameters
fullDefault=True

# Gaussian parameters
title = "title"
charge = "0"
multiplicity = "1"
#functional = "mpw1pw91"
functional = "mpw1pw91"
#basis = "6-31g(d)"
basis = "6-31g(d)"
#nprocshared = "8"
nprocshared=""
#mem="16GB"
mem=""
#addnl_route = "scrf=(solvent=generic,read)"
addnl_route=''

# Additional gaussian lines (at end of file after xyz coordinate setion)
#extra_lines = "Eps=33 \nEpsinf=2.15 \n" #Put \n characters as new lines
extra_lines=''

jt=0 # 0 for both, 1 for nics, 2 for mcbo 

# ∫NICSzz parameters
tot_points = 21
bq_max = 3
bq_min = -3
bq_spacing = 0.1

bq_points = np.arange(bq_min, bq_max+0.00001, bq_spacing)

# Input file name options
directory_name = '_from_'
outfilename = '_REPLACE.gjf'

# -------------------------------------------------------------------------- #
# END OF PARAMETER INITIALIZATION

def init_vars():
    
    return (infilename, fullDefault, title, charge, multiplicity, functional, 
            basis, nprocshared, mem, addnl_route, extra_lines, jt, tot_points, bq_max, bq_min, bq_spacing,
            bq_points, directory_name, outfilename)

def get_lines(infilename):
    lines = None
    with open(infilename) as f:
        lines = f.readlines()
    return lines
    

def get_coords(lines, lt):
    flip = None
    irc_coords = []
    for i in range(len(lines)):
        line = lines[i]

        if "Input orientation:" in line:  # Finding indicator of coordinates in log file
            curr_i = i + 5  # Coordinates begin 5 lines after this keyword
            coords = []
            if "Total number" in lines[i-2]:
                continue
            else:
                while '------' not in lines[curr_i]:  # -- indicates end of coordinates
                    rm_space = " ".join(lines[curr_i].split())
                    stripped = np.asarray(rm_space.split(' '))
                    curr_coords = stripped[[-5, -3, -2, -1]]
                    coords.append(curr_coords)
                    curr_i += 1
                    
                irc_coords.append(coords)

        elif 'ning calculation of the REVERSE path' in line:  # Find where IRC flips direction
            flip = len(irc_coords)
            
    irc_coords = np.asarray(irc_coords).astype(float)
    
    if lt == 0: # IRC

        ts_step = len(irc_coords) - flip
        organized_irc_coords = np.concatenate((np.flipud(irc_coords[flip:]), irc_coords[:flip]))
        print("TS_step on IRC is: ", ts_step)
        return organized_irc_coords, ts_step
        return irc_coords, 1
    
    elif lt == 1: # TRAJ
        return irc_coords, 1
    
    elif lt == 2: # OPT
        return irc_coords, 1
        
def get_log_type(lines): #Returns: 0=IRC, 1=TRAJ, 2=OPT
    for line in lines:
        if "IRC-IRC-" in line:
            print("Log file is an IRC.")
            return 0
        
        elif "TRJ-TRJ-" in line:
            print("Log file is a trajectory.")
            return 1
        
    print("Log file not IRC or TRAJ, assuming is an optimization or singlepoint")
    return 2
        
def get_func_basis(lst):
    # Iterate through the list to find an element containing "/"
    for element in lst:
        if "/" in element:
            print(element)
            func_basis = element.split("/")
            return func_basis
    
    # If no element contains "/", prompt the user for data
    func_basis = input("Func/Basis not found in route line. Please provide input with 'func/basis': ").split("/")
    return func_basis

def get_func_bas_chrg_mult(lines):
    f=None
    b=None
    c=None
    m=None
    route_line=""
    for i in range(len(lines)):
        line = lines[i]
        if " Charge = " in line:
            splitline = ' '.join(line.split()).split(' ')
            c = splitline[2]
            m = splitline[-1]
            print(f'Charge and Multiplicty found. Charge: {c}, Multiplicity: {m}')
            
        elif " #" in line:
            curr_line = lines[i]
            while "---" not in curr_line:
                route_line += curr_line.strip()
                i += 1
                curr_line = lines[i]
        if route_line != "" and all(var is not None for var in (c, m)):
            print(f'Route line found: {route_line}')
            break
    splitline = route_line.split(' ')
    f,b=get_func_basis(splitline)              


    if None not in (f,b,c,m): # If each variable has been given a value
        print("Func, Basis, Charge, Mult: {} {} {} {}".format(f,b,c,m))
        return (f,b,c,m)


def prompt_jt( jt): # Returns the jobtype to create, none=both , 0=cancel, 1=nics, 2=mcbo
    
    try:
        jt = int(input('Type 0 to cancel, 1 to create ∫NICSzz input files, 2 to create MCBO input files: '))
    except:
        print('Please input a valid option')
    
    if jt=="":
        print("Creating input files for NICS and MCBO analysis...")
        
    elif jt=="0":
        print("Canceling...")
        sys.exit()
        
    elif jt=="1":
        print("jt = 1, creating ∫NICSzz input files...")
        
    elif jt=="2":
        print("jt = 2, creating sp jobs with checkpoint files for Multiwfn analysis...")
        
    else:
        print("Invalid job option, try again...")
        sys.exit()
    return jt

def set_out_dir_name(jt, bq_min, bq_max, spacing):
    if jt=="":
        return "temp"
    

def prompt_bq_spacing(bq_min, bq_max, spacing, bq_points):
    print('Bq atoms to sample NICS are placed as follows:')
    print('Above ring: {}Å'.format(bq_max))
    print('Below ring: {}Å,'.format(bq_min))
    print('Spacing: {}Å'.format(bq_spacing))
    print('Bq points are {}'.format(bq_points))
    uin = input("Press enter to proceed with default, 0 to cancel, 1 to do manual input. (Parameters can be adjusted in executable file):")
    
    if uin == "":
        print("Proceeding with default options...")
        
    elif uin=="0":
        print("Canceling...")
        sys.exit()
        
    elif uin=="1":
        uin_tot_points = input("Input files by irc step, press enter to use default points, or enter step # of irc you want a NICSzz analysis of (i.e 1,5,10,15,20)")
        if uin_tot_points == "":
            print("Using default irc points...")
        else:
            try:
                split_points = uin_tot_points.split(',')
                bq_points = [int(x) for x in split_points]
                print("Manually enetered irc steps to analyze:")
                print(bq_points)
                return bq_points
                
            except:
                print("Invalid value, exiting...")
                e = sys.exc_info()[0]
                print("ERROR:")
                print(e)
                sys.exit()
    return bq_points
    

def prompt_rxn_points(rxn_points):
    uin = input("Input files by irc step, press enter to use default points, or enter step # of irc you want a NICSzz analysis of (i.e 1,5,10,15,20)")
    if uin == "":
        print("Using default irc points...")
        return rxn_points
    else:
        try:
            split_points = uin.split(',')
            points = [int(x) for x in split_points]
            print("Manually enetered irc steps to analyze:")
            print(points)
            return(points)
            
        except:
            print("Invalid value, exiting...")
            sys.exit()

def prompt_arom_atoms(): 
    # Requests user for aromatic atom indices.
    # Two fragments must be specified for a Diels-Alder.
    # Otherwise SVD will fail to define a good normal vector
                         
    arom_atoms = input("Please type the indicies of the aromatic carbons (ie. 1,2,3,4...):\
                      If diels-alder type reaction, enter the indices of aromatic atoms of two fragments seperated by a space (i.e 1,2,3,4 5,6): ")
    open("arom_atoms.txt","w").write(arom_atoms)
        
    
    frag1 = None
    frag2 = None                                                                                                                            
    if ' ' in arom_atoms:
        frag1, frag2 = arom_atoms.split(' ')
        frag1 = np.asarray([int(x) for x in frag1.split(',')])
        frag2 = np.asarray([int(x) for x in frag2.split(',')])
    else:
        frag1 = np.asarray([int(x) for x in arom_atoms.split(',')])
    return frag1, frag2
    
def split_irc_steps(tot_points, ts_step, max_step, lt):
    if lt == 0: # IRC > split to front and back of TS
        left_irc = np.linspace(1, ts_step, int(np.ceil(tot_points / 2)), dtype=int)
        print(f'running: np.linspace({ts_step},{max_step},{int(np.ceil(tot_points/2))},dtype=int')
        right_irc = np.linspace(ts_step, max_step, int(np.ceil(tot_points / 2)), dtype=int)
        rxn_points = np.concatenate((left_irc[:-1], right_irc))
        return rxn_points
    
    elif lt == 1: # TRAJ, no need to split
        rxn_points = np.linspace(1, max_step, tot_points, dtype=int)
        return rxn_points
    
    elif lt ==2: # OPT OR SP
        rxn_points=np.asarray([-1])
        print(f'rxn_points:{rxn_points}')
        return rxn_points

def set_dir_name(dir_name, jt):
    if jt == 0:
        dir_name='nics_mcbo_from_irc'
    if jt == 1:
        dir_name='nics_from_irc'
    if jt == 2:
        dir_name='mcbo_from_irc'
    #return dir_name
    return "peri_irc"

def create_gaussian_job(bq_points, rxn_points, irc_coords, jt, directory_name, arom_atoms, nprocshared, mem, functional, basis, title, charge, multiplicity, addnl_route, extra_lines):
        
    try:  # Tries creating new directory to place .com files in
        os.mkdir(directory_name)  # Making directory for files
    except FileExistsError:
        print('directory: ', directory_name, ' already exists')
        
    for n in rxn_points: # Loop through each point being sampled in IRC/TRAJ
        coords = irc_coords[n-1][:][:] # -1 index bc indexing starts at 0 instead of 1
        nics_coords = coords.tolist() # Make a list where Bq atoms will be added
        for nc in nics_coords:
            nc = np.asarray(nc)
        #print(nics_coords)

        for i in range(len(nics_coords)): # Loop through existing coordinates
            nics_coords[i][0] = pt.elements[nics_coords[i][0]]
            
            if arom_atoms[1] is None: #If only one fragment was specified
                xyz = coords[arom_atoms[0] - 1][:, 1:]  # Just xyz coords of aromatic carbons
                cent = centroid(xyz)  # centroid xyz
                nvec = norm_vect(xyz)  # vector normal to best fit plane to aromatic carbons
                
            else: # If two fragments were specified:
                xyz1 = coords[arom_atoms[0] - 1][:, 1:]
                xyz2 = coords[arom_atoms[1] - 1][:, 1:]
                combined_frags = np.concatenate((arom_atoms[0], arom_atoms[1]), axis=None)
                xyz = coords[combined_frags - 1][:, 1:]
                
                cent1 = centroid(xyz1)
                cent = centroid(xyz)
                
                #Find normal vector for fragment1
                nvec1 = norm_vect(xyz1)
                
                #Shift xyz2 so coplanar with xyz1
                xyz2_shifted=[]
                for coord in xyz2:
                    dist = cent1-coord
                    xyz2_shift = project_vec(dist, nvec1)
                    xyz2_shifted.append(xyz2_shift+coord)
                nvec = norm_vect(np.vstack([xyz1,xyz2_shifted]))
        
        for nics_index in bq_points:
            nc=cent - nics_index * nvec
            nics_coords.append(['Bq', str(nc)])
            
        
        output_template = get_output_template(n, nprocshared, mem, jt, functional, basis, title, charge, multiplicity, addnl_route)
        if n==-1:
            n=1
        outfilename=directory_name+'_REPLACE.gjf'
        with open(directory_name + '/' + outfilename.replace('REPLACE', str(n).zfill(4)), 'w') as f:
            for line in output_template:
                f.write(line)
            for line in nics_coords:
                new_line = ("%s\n" % line).replace('[', '').replace(']', '').replace('[', '').replace(',', '').replace('\'', '').replace('\n', '')
                l = re.sub(' +', ' ', new_line)
                #print(l)
                spline = l.split(' ')
                fline = ''
                #print(spline)
                for i, x in enumerate(spline):
                    if x == '':
                        pass
                    elif i != 0:
                        y=float(x)
                        fline += f'{y:.6f}'
                    else:
                        fline += x
                    fline += ' '
                fline+='\n'
                f.write(fline)
                
                #f.write(new_line)
            f.write('\n')
            f.write(extra_lines)
            f.write('\n')
        
def get_coords_from_line(line):
    """
    Parses out atomic number and xyz coordinates from line
    :param line: string
    :return: 1x4 array, [atomic number, x, y, z]
    """
    rm_space = " ".join(line.split())
    stripped = np.asarray(rm_space.split(' '))
    fstripped=[]
    for x in stripped[-5,-3,-2,-1]:
        fstripped.append(float(x))
    return np.asarray(fstripped)
    
#----------------------------------------------------------------------------#
# MATH FUNCITONS for defining normal vector to aromatic atoms
#----------------------------------------------------------------------------#

def project_vec(v1,v2):
    return np.dot(v1, v2)/np.linalg.norm(v2)**2*v2

def centroid(points):
    """
    Calculates the centroid of a set of xyz points
    :param points: Nx3 2D Array of xyz coordinates of N atoms
    :return: 1x3 array, [x,y,z], coordinates of centroid
    """
    length = len(points)
    sum_x = np.sum(points[:, 0])
    sum_y = np.sum(points[:, 1])
    sum_z = np.sum(points[:, 2])
    return [sum_x / length, sum_y / length, sum_z / length]


def norm_vect(points):
    """
    Finds vector that is normal to the plane that best fits
        the points passed in. Finds through SVD.
    :param points: Nx3 2D Array of xyz coordinates of N atoms
    :return: 1x3 array, [x,y,z], vector normal to the plane
        that best fits the points
    """
    points_m_centroid = points - centroid(points)
    u, s, vh = np.linalg.svd(points_m_centroid)
    v=vh.transpose()
    nvec = v[:, 2]
    #print('normal vector')
    #print(nvec)
    if nvec[2] < 0:  # Ensure normal vector is always in +z direction
        nvec *= -1
    return nvec

def get_output_template(rxn_step, nprocshared, mem, jt, functional, basis, title, charge, multiplicity, addnl_route):
    outfilename = ""
    output_template = ""
    chkline = ""
    if rxn_step==-1:
        rxn_step=1
        print("reset rxn step")
    if jt == 0:
        jobtype = "nmr"
        outfilename="peri_irc_REPLACE.com"
        chkline = '%chk=' + outfilename.replace('REPLACE', (str(rxn_step).zfill(4)))[:-3] + 'chk \n'
    if jt == 1:
        jobtype = "nmr"
        outfilename="peri_irc_REPLACE.com"
        
    elif jt == 2:
        jobtype =""
        outfilename="peri_irc_REPLACE.com"
        chkline = '%chk=' + outfilename.replace('REPLACE', (str(rxn_step).zfill(4)))[:-3] + 'chk \n'
    
    if nprocshared != "":
        output_template += '%NProcShared=' + nprocshared + '\n'
    
    if mem != "":
        output_template += '%Mem=' + mem + '\n'
    
    if jt==2 or jt==0: # Insert chkpoint line for MCBO analysis
        output_template += chkline
        
    # Constructucting gaussian input template
    gau_essentials = '# ' + jobtype + ' ' + functional + '/' + basis + ' ' + addnl_route + '\n', \
                       '\n', \
                       title + '\n', \
                       '\n', \
                       charge + ' ' + multiplicity, \
                       '\n'
                       
    gau_essentials = ''.join(gau_essentials)
    
    output_template += gau_essentials
        
    return output_template

def plot_atoms(coords):
    print('IN PLOT FUNCTION')
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    for i in range(len(coords)):
        atom=coords[i]
        atom_name = atom[0]
        atom_xyz = atom[1:]
        
        i=0
        if atom_name == 6.:
            c = "gray"
        elif atom_name == 1.:
            c = "white"
        elif atom_name == 8.:
            c = "red"
        elif atom_name == 7.:
            c = "blue"
        else:
            c = "purple"
            i+=1
            
        ax.scatter(*atom_xyz, color = c)
        ax.text(*atom_xyz, i, color = c)
        ax.scatter(*atom_xyz,s=atom_name*50, color = c)
    print("FINISHED PLOT FUNCTION")
        
#----------------------------------------------------------------------------#
# BELOW FUNCTION ONLY COMPATIBLE WITH IPYTHON, USE plot_atoms FUNCTION INSTEAD
# IF IPYTHON IS NOT AVAILABLE
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#    
def iplot_atoms(coords):
    plt.ion
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    for atom in coords:
        c=""
        atom_name = atom[0]
        atom_xyz = atom[1:]
        if atom_name == 6.:
            c = "gray"
        elif atom_name == 1.:
            c = "white"
        elif atom_name == 8.:
            c = "red"
        elif atom_name == 7.:
            c = "blue"
        else:
            c = "purple"
        ax.scatter(*atom_xyz, color = c)
        ax.text(*atom_xyz, atom_name, color = c)
        ax.scatter(*atom_xyz,s=atom_name*50, color = c)

#----------------------------------------------------------------------------#
# Runs perixn job preparation
#----------------------------------------------------------------------------#
def run(infilename=None):
    np.set_printoptions(suppress=True)
    #Initializing preset variables
    infilename_def, fullDefault, title, c, m, f, b, nprocshared, mem, addnl_route, extra_lines, jt, tot_points, bq_max, \
    bq_min, bq_spacing, bq_points, directory_name, outfilename = init_vars()
    
    if infilename is None: # If no file is specified, use default infilename
        infilename=infilename_def
    
    print("infilename:", infilename)
    lines = get_lines(infilename) # Get lines of log file
    
    f,b,c,m = get_func_bas_chrg_mult(lines) # Get functional, basis set, charge
                                            # and mulitplicity from log file
    lt = get_log_type(lines) # Get type of log file, IRC, TRAJ, OPT
    
    coords, ts_step = get_coords(lines, lt) # Get coordinates of IRC/TRAJ/OPT
                                            # Also returns the step # of ts
    #plot_atoms(coords[5])
    print(coords[:,0])
    print(len(coords))
    bq_points = np.arange(bq_min, bq_max+0.00001, bq_spacing)
    #rxn_points=[32,33,34,36,37,38]
    
    rxn_points = split_irc_steps(tot_points, ts_step, len(coords), lt)
    print(f'rxn_points: {rxn_points}')
    
    if not fullDefault: # Set fullDefault=False to allow user override/check
        # Specify mNICS, MCBO or both job preparation
        jt = prompt_jt()
        
        # Array of bq atoms to sample nics data
        bq_points = prompt_bq_spacing(bq_min, bq_max, bq_spacing)
        
        # Array of points to sample from IRC/TRAJ
        rxn_points = prompt_rxn_points()
    
    arom_atoms = prompt_arom_atoms()
    
    directory_name = set_dir_name(directory_name, jt)
    
    create_gaussian_job(bq_points, rxn_points, coords, jt, directory_name, arom_atoms, nprocshared, mem, f, b, title, c, m, addnl_route, extra_lines)
    
    shutil.copyfile("./"+infilename, "./peri_irc/"+infilename)

def main(argv):
    if not argv:
        filename=input("Please type the name of the log file.\n")
        run(filename)
    else:
        run(argv[0])
        
    
if __name__ == "__main__":
    main(sys.argv[1:])

