#!/usr/bin/python

"""
This script takes Gaussian output and uses the provided params file as a 'template'
to put in better internal coordinates.
"""

if len(sys.argv) < 4:
    print("Run this script like so:")
    print("python internal_coordinates.py $mol.out $ref.pdb $params")
    print("Doing so will take the Gaussian output at $mol.out, use the pdb file in the ")
    print("present working directory at $ref.pdb to map them to atom names, and use the ")
    print("params file $params as a template to put in updated ")
    print("coordinates.")
    print("")
    print("In the future, with any luck, we'll have a way to impose RESP charge fitting")
    print("results here too...")


ROSETTA = os.environ['ROSETTA']
# Set your Rosetta directory here, if for some reason you can't use environment variables.
#ROSETTA = "/path/to/Rosetta/"

import sys
import numpy as np
import numpy.linalg as la
import os

# MATH
def dist( p1, p2 ):
    return ( ( p1[0]-p2[0] )**2.0 + ( p1[1]-p2[1] )**2.0 + ( p1[2]-p2[2] )**2.0 ) ** (0.5)

def ang( p1, p2, p3 ):
    v1 = p1-p2 
    v2 = p3-p2 
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return 180.0-np.degrees(np.arctan2(sinang, cosang))

def tors( p1, p2, p3, p4 ):
    b = np.array([p2-p1, p3-p2, p4-p3])
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return -1.0 * np.degrees(np.arctan2( y, x ))

# PDB PROCESSING
def atomno( line ):
    return int(line[6:12].strip())

def atomname( line ):
    return line[12:16].strip()


# START
def get_mulliken_charges( fn ):
    charge_set = []
    f = open(fn)
    lines = f.readlines()
    f.close()
    
    for i in xrange( len(lines) ):
        line = lines[i]
        if line[0:len(" Mulliken charges:")] == " Mulliken charges:":
            # next natoms lines
            chg = []
            while line[0:3] != "Sum":
                i += 1
                line = lines[i].strip()
                if len( line.split() ) < 3: continue
                if line[0:3] == "Sum": break
                print(line)
                chg.append( float(line.split()[2]) )
            charge_set.append(chg)

    return charge_set

def process_cmd_line( argv ):
    prefix = ""
    gauss_fn = sys.argv[1]
    pdb_fn = sys.argv[2]

    import cclib
    from cclib.parser import Gaussian
    #from cclib.method import CSPA #MPA
    myfile = Gaussian(gauss_fn)
    data = myfile.parse()
    opt_coords = data.atomcoords[-1] # last coords from the optimization
    
    charges = get_mulliken_charges( fn )[-1]
    #analysis = CSPA(myfile)#MPA(myfile)
    #if analysis.calculate():
    #    charges = analysis.fragcharges # last coords from the optimization
    #    print(charges)
    
    paramslines = [] 
    try:
        with open(sys.argv[3]) as f:
            paramslines = f.readlines()
    except:
        with open(ROSETTA + "/main/database/chemical/residue_type_sets/fa_standard/residue_types/{}.params".format(params)) as f:
            paramslines = f.readlines()
    
    pdb = []
    with open(pdb_fn) as f:
        pdb = f.readlines()
    
    return opt_coords, charges, pdb, paramslines


opt_coords, charges, pdb, params = process_cmd_line( sys.argv )

def point_from_line( line ):
    objs = line.split()
    return Point(float(objs[0]), float(objs[1]), float(objs[2]) )

name_to_num = { atomname(line): atomno(line) for line in pdb if line[0:5] == "ATOM " }
num_to_name = { atomno(line): atomname(line) for line in pdb if line[0:5] == "ATOM " }
print(num_to_name)
for num in num_to_name.keys():
    if num_to_name[num] == "YP" or num_to_name[num] == "NM":
        num_to_name[num] = "UPPER"
    elif num_to_name[num] == "O3P" or num_to_name[num] == "CO":
        num_to_name[num] = "LOWER"

coords = { num_to_name[num] : opt_coords[num-1] for num in num_to_name.keys() }
print(coords)
# Read each line in params. If it contains ICOORs, update appropriately.
# Otherwise, leave it alone.
counter = 0

out_params = open("out.params", "w")
for line in params:
    #print(line)
    if line[0:5] != "ICOOR" and line[0:5] != "ATOM ":
        out_params.write(line)
        continue
    #print(line    )
    if line[0:5] == "ATOM ":
        objs = line.strip().split()
        if objs[2] == "VIRT": out_params.write( line ) 
        else: out_params.write("{} {:>4s} {:>4s} {:7.6f}\n".format(line[0:9], objs[2], objs[3], charges[name_to_num[objs[1]]-1]))
        continue
    #print(line)
    objs = line.split()

    counter = counter + 1
    name = objs[1]

    # If name isn't in coords--say, if it's a virt like HO2' in O-Me NTs--skip it.
    if not name in coords: 
        out_params.write(line)
        continue

    stub1 = objs[5]
    stub2 = objs[6]
    stub3 = objs[7]

    phi = 0
    theta = 0
    d = 0

    if counter >= 2:
        d = dist( coords[name], coords[stub1] )
    if counter >= 3:
        theta = ang( coords[name], coords[stub1], coords[stub2] )
    if counter >= 4:
        phi = tors( coords[name], coords[stub1], coords[stub2], coords[stub3] )
    
    out_params.write("ICOOR_INTERNAL    %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (name, phi, theta, d, stub1, stub2, stub3 ) )
out_params.close()
