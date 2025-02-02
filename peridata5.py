#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 11:34:42 2021

@author: willdesnoo
"""

import numpy as np
import sys
import os
import xlsxwriter as xw
import subprocess
#import periparser


hartree2kcal=627.509608
script_path = os.path.dirname(os.path.realpath(__file__))
multiwfn_path = '/Applications/ws/Multiwfn'
mcbo_script_path = script_path + '/nmcbo_irc_v4'
elf_script_path = script_path + '/elf_irc.sh'
formchk_all_script_path = script_path + '/formchk_all'

def find_irc_file(path):
    lfs=[]
    print("No original log file provided, searching current directory...")
    for file in os.listdir(path):
        if file.endswith(".log"):
            lfs.append(file)
            
    print("Found log files:")
    print("----------------")
    print("INDEX ---- FILENAME")
    for i in range(len(lfs)):
        print("  {}   ---- {}".format(i, lfs[i]))
        
    uin=int(input("Please type the index of the original log file of which the nics"
              " and nMCBO were derived from...\n"))
    
    print("Chosen original log file is: {}".format(lfs[uin]))
    
    return(lfs[uin])

def get_files_in_dir(cwd,formchk_all_script_path, verbose=False):
    path=cwd
    log_files = []
    fchk_files = []
    
    for file in os.listdir(path):
        if file.endswith(".log"):
            log_files.append(file)
        elif file.endswith(".fchk"):
            fchk_files.append(file)
    if verbose:        
        print("For mNICS: Found {} log files...".format(len(log_files)))
        print("For nMCBO: Found {} fchk files...".format(len(fchk_files)))
    
    if len(fchk_files)==0:
        print("No .fchk files detected, running formchk_all to create .fchk files.")
        subprocess.check_output([formchk_all_script_path])
        
        
    return log_files, fchk_files

def get_nics(lfs):
    indices=[]
    irc_step=[]
    full_bq_list=[]
    # Find all log files
    for f in lfs:
        bqlist=[]
        try:
            irc_step.append(int(f.split('.')[0].split('_')[-1])) # Get irc step from name of input file
        except(ValueError):
            print("File name does not include IRC step, continuing without IRC step labeling")
        indices, bqlist = find_Bq_zz(f) # Find Bq atoms
        full_bq_list.append(bqlist)
    full_bq_array=np.asarray(full_bq_list)
    #full_bq_array = np.asarray([x for _, x in sorted(zip(irc_step, full_bq_array))])
    return irc_step, indices, full_bq_array

def get_mcbo(mcbo_scripts_path,verbose=False):
    nmcbo_raw=None
    if os.path.isfile("arom_atoms.txt"):
        with open("arom_atoms.txt","r") as file: arg=','.join(file.readlines())
        print(arg)
    else:
        arg = input("Input indicies of aromatic atoms for nMCBO analysis (i.e 1,2,3,4,5,6).\n")
        with open("arom_atoms.txt","w") as file: file.write(arg)
    if verbose:
        print("Running MCBO command as: {} {}".format(mcbo_script_path, arg))
    try:
        nmcbo_raw = subprocess.check_output([mcbo_script_path, arg])
    except subprocess.CalledProcessError as exc:
        nmcbo = exc.output.decode("utf-8")
        return nmcbo
    else:
        nmcbo = nmcbo_raw.decode("utf-8")
        if verbose:
            print(nmcbo)
        return nmcbo.split('\n')[:-1]
    '''
    except Exception as err:
        print("Error when calculating nMCBO values...")
        print(err)
        print("Try manually using mcbo script to get MCBO values...")
        sys.exit()
    '''
    

def get_elf(filename):
    try:
        elf_out=subprocess.check_output([elf_script_path, filename])
        return elf_out
    
    except Exception as err:
        print("Error when calculating ELF values at critical points...")
        print(err)
        print("Try manually using ELF script to get ELF values.")
    
def find_Bq_zz(filename):
    bq_list = []
    indices = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if 'Bq   Isotropic' in line:
                index = lines[i].split(' ')[5]
                zz = lines[i + 3].split(' ')[-1]
                bq_list.append(-float(zz))
                indices.append(index)
        f.close()
    return indices, np.asarray(bq_list)

def save_nics_data(workbook, bq_array, irc_step, bq_placement=np.arange(-3, 3.0001, 0.1)):
    dim1 = len(bq_array)
    dim2 = len(bq_array[0])
    
    ws1 = workbook.add_worksheet("mNICSzz")
    ws1.write(0, 0, 'Reaction Step')

    ws1.write(0, 1, "mNICSzz")

    for i in range(dim2):
        ws1.write(0, i + 2, 'NICS(%0.1f)zz' % bq_placement[i])
        
    # Writing row names (irc step)
    for i in range(dim1):
        if len(irc_step) != 0:
            ws1.write(i+1, 0, irc_step[i])
        else:
            ws1.write(i+1, 0 , 'curr_structure')

    num_format = workbook.add_format({'num_format': '0.00'})
    for i in range(dim1):
        curr_nics_array = bq_array[i]
        for j in range(dim2):
            ws1.write(i+1, 1, float(curr_nics_array.mean()), num_format)
            ws1.write(i+1, j + 2, float(curr_nics_array[j]), num_format)
    

def get_energies(filename):
    print("Reading coordinates from irc...")
    flip=False
    fpath=[]
    rpath=[]
    energy=0
    ts=None
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            #if 'Energy From Chk' in line or 'Energies reported relative to the TS energy of' in line:
            if 'Energies reported relative to the TS energy of' in line:
                ts = [0, 0, float(line.split()[-1])]
                rpath.append(ts)
            if 'SCF Done:' in line:
                energy = float(line.split()[4])
                
            if "Pt" in line:  # Finding indicator of coordinates in log file
                split_line=lines[i+1].split()
                step = int(split_line[2])
                
                split_line=lines[i+3].split()
                irc_step=float(split_line[-1])
                
                if flip:
                    if int(-1*step) in column(rpath,0):
                        print('duplicate!')
                    else:
                        rpath.append([-1*step, -1*irc_step, energy])
                else:
                    fpath.append([step, irc_step, energy])
    
            elif 'Begin' in line:  # Find where IRC flips direction
                flip=True
    fullpath = np.asarray(rpath+fpath)
    try:
        irc1=min(fullpath[:,0])
        fullpath[:,0] -= (irc1-1)
    except Exception:
        pass
    #return np.asarray(fullpath)
    if len(fullpath) > 0:
        energies=fullpath[np.argsort(fullpath[:,0])]
        return energies
    else:
        return [energy]

def save_energies(ws, energies):
    ws.write(0,0,"Reaction Step")
    ws.write(0,1, "Intrinsic Reaction Coordinate")
    ws.write(0,2, "Electronic Energy (Hartree)")
    ws.write(0,3, "Relative Electronic Energy (kcal/mol)")    
    
    if len(energies) > 1:
        e1=energies[0][2]
        for i in range(len(energies)):
            curr_e = energies[i]
            for j in range(len(curr_e)):
                ws.write(i+1, j, curr_e[j])
            ws.write(i+1, len(curr_e+1), hartree2kcal*(curr_e[2]- e1))
    else:
        e1=float(energies[0])
        ws.write(1,1,e1)

def save_mcbo0(ws, mcbo):
    ws.write(0,5,"nMCBO")
    for i,bo in enumerate(mcbo):
        ws.write(i+1, 5, mcbo[i])

def save_mnics0(ws, mnics):
    ws.write(0,4,"mnics")
    for i,mn in enumerate(mnics):
        ws.write(i+1,4,mnics[i])
        
def save_nics_with_energies(ws,wb,nics,mnics0,bq_placement=np.arange(-3,3.001,0.1),verbose=True):
    wid=len(nics[0])

    for i in range(wid):
        ws.write(0, i + 6, 'NICS(%0.1f)zz' % bq_placement[i])
    
    num_format = wb.add_format({'num_format': '0.00'})
    c=0
    if verbose:
        print("NICS")
        print(np.shape(nics))
        print(f"len(mnics0): {len(mnics0)}")
        print(f'nicsShape: {np.shape(nics)}')
        print(mnics0)
        print(wid)
    for i in range(len(mnics0)):
        try:
            curr_nics_array = nics[c]
        except Exception as e:
            print('caught exception')
            print(e)
            print(i)
        if mnics0[i] != None:
            for j in range(wid):
                ws.write(i+1, j + 6, float(curr_nics_array[j]), num_format)
            c+=1
            
            
def save_dictionary(workbook, d):
    worksheet = workbook.add_worksheet("The Kitchen Sink")
    row = 0
    col = 0
    
    order=sorted(d.keys())
    for key in order:
        row += 1
        worksheet.write(row, col,     key)
        i=1
        for item in d[key]:
            for it in item:
                worksheet.write(row, col + i, item)
                i += 1
    
def save_elf(wb, elfs):
    print('Saving elfs...')
    ws = wb.add_worksheet()
    for i,elf in enumerate(elfs):
        row = 0
        col = i+2
        rowi = 0
        for val in elf:
            sval = val.split(',')
            if i == 0:
                ws.write(rowi, 0, sval[0])
                ws.write(rowi, 1, str(sval[1]+','+sval[2]))
                rowi+=1
            #print(f'val: {val}')
            ws.write(row, col, float(sval[-1]))
            row+=1
        
def save_bond_elf(wb,elfs):
    print('Saving elfs...')
    for i,elf in enumerate(elfs):
        ws = wb.add_worksheet()
        row = 0
        for val in elf:
            sval = val.split(',')
            if ("-1" in sval[2]): #If bond CP
                #print(sval)
                ws.write(row, 0, sval[0])
                ws.write(row, 1, str(sval[1]+','+sval[2]))
                ws.write(row, 2, float(sval[3]))
                row+=1

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def column(matrix, i):
    return [row[i] for row in matrix]
                
def run(ircfile,plot=True):
    workbook = xw.Workbook('nics_mcbo_data.xlsx')
    ws1=workbook.add_worksheet('IRC Energies')
    lfs, ffs = get_files_in_dir(os.getcwd(),formchk_all_script_path)
    
    print(f'IRC Log File: {ircfile}')
    if ircfile in lfs:
        lfs.remove(ircfile)
    
    full_bq_array=[]
    energies = get_energies(ircfile)
    
    irc_step, indices, full_bq_array=get_nics(lfs)
    irc_step, full_bq_array = zip(*sorted(zip(irc_step, full_bq_array)))
    full_bq_array=np.asarray(full_bq_array)
    #mnics=np.average(full_bq_array[:,21:41], axis=1)
    mnics=np.average(full_bq_array, axis=1)
    save_nics_data(workbook, np.asarray(full_bq_array), irc_step)
    mcbo=get_mcbo(mcbo_script_path)

    full_irc_step=np.arange(1,len(energies)+1)
    
    mcbo0=[]
    mnics0=[]
    c=0
    for i in range(len(full_irc_step)):
        if (i+1) in irc_step:
            mcbo0.append(float(mcbo[c]))
            mnics0.append(mnics[c])
            c+=1
        else:
            mcbo0.append(None)
            mnics0.append(None)
    save_mcbo0(ws1, mcbo0)
    save_mnics0(ws1, mnics0)
    save_energies(ws1,energies)
    save_nics_with_energies(ws1, workbook, np.asarray(full_bq_array),mnics0)
    
    workbook.close()
    

    if plot:
        import periplot as pp
        energiesH=np.asarray(column(energies, 2))
        firstPt=np.asarray(column(energies, 2))[0]
        energiesKcalMol=[hartree2kcal*(x-firstPt) for x in energiesH]

        pp.plot_data(np.asarray(column(energies, 0)), energiesKcalMol, np.asarray(mnics0,dtype=float), np.asarray(mcbo0), np.asarray((column(energies, 1))))
        pp.plotNICSHM(np.asarray(column(energies, 0)), energiesKcalMol,np.asarray(full_bq_array),irc_step)
        
def write_dat(data,out_name='peri_data.out'):
    with open(out_name) as file:
        file.write(data)

def main(argv):
    # Initiate arguement parser for ash-like options
    
    if not argv:
        argv = find_irc_file(os.getcwd())
    run(argv)
    #test()
        
    
if __name__ == "__main__":
    main(sys.argv[1:])

