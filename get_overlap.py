"""
ADF Input File Generator
Author: Dennis Svatunek


Description:
This script reads orbital overlaps from EDA calculations
"""

import os
import csv
import numpy as np
from scm.plams import *
init()


def GetAtom(kf):
   # preparations: get geometry info, atomic symbols in the right order
   nAtoms = kf.read('Geometry', 'nr of atoms')
   aAtoms = kf.read('Geometry', 'fragment and atomtype index')[nAtoms:]
   xAtoms = str(kf.read('Geometry', 'atomtype')).split()
   oAtoms = kf.read('Geometry', 'atom order index')[nAtoms:]
   sAtoms = [xAtoms[order-1] for numb, order in sorted(zip(oAtoms, aAtoms))]
   
   xyz = kf.read('Geometry', 'xyz')
   # Convert the list to a numpy array before reshaping
   xyz_array = np.array(xyz)
   coordinates = xyz_array.reshape((nAtoms, 3))
   
   # Convert coordinates from Bohr to Angstrom
   coordinates = coordinates * 0.52917721092

   
   return nAtoms, sAtoms, coordinates

def calculate_distance(coord1, coord2):
    """
    Calculate the distance between two points in 3D space.
    
    Parameters
    ----------
    coord1, coord2 : numpy.ndarray
        The x, y, z coordinates of each point.
    
    Returns
    -------
    float
        The distance between the points.
    """
    return np.linalg.norm(coord1 - coord2)

def calculate_angle(coord1, coord2, coord3):
    """
    Calculate the angle (in degrees) between three points in 3D space.
    
    Parameters
    ----------
    coord1, coord2, coord3 : numpy.ndarray
        The x, y, z coordinates of each point.
    
    Returns
    -------
    float
        The angle in degrees.
    """
    v1 = coord1 - coord2
    v2 = coord3 - coord2
    cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)
    
    
def GetFaIrrep(complexResult):
    irrepType = str(complexResult.read('Symmetry', 'symlab')).split()
    irrepOrbNum = complexResult.read('Symmetry', 'norb')

    if len(irrepType) == 1:
        irreporbNum = [irrepOrbNum]
    else:
        irreporbNum = irrepOrbNum

    faIrrepone = [[irrep for _ in range(number)] for irrep, number in zip(irrepType, irreporbNum)]
    return [irrep for sublist in faIrrepone for irrep in sublist]

def ReadOverlap(complexResult, faIrrep, index_1, index_2):
    faOrb = complexResult.read('SFOs', 'isfo')
    maxIndex = max(faOrb[index_1], faOrb[index_2])
    minIndex = min(faOrb[index_1], faOrb[index_2])
    index = maxIndex * (maxIndex - 1) // 2 + minIndex - 1

    if faIrrep[index_1] == faIrrep[index_2]:
        overlap_matrix = complexResult.read(faIrrep[index_1], 'S-CoreSFO')
        #print(overlap_matrix)
        return abs(overlap_matrix[index])
    else:
        return 0
    
def GetOrbitalIndex(complexResult, Frag, MO):
      orbIndex = 0
      orbEnergy=complexResult.read('SFOs', 'energy')
      orbFragment = complexResult.read('SFOs', 'fragment')
      fragIrrep = str(complexResult.read('SFOs', 'subspecies')).split()
      fragOrb  = complexResult.read('SFOs', 'ifo')
      for i in range(len(orbEnergy)):
        if orbFragment[i] == Frag  and  fragIrrep[i] == 'A' and fragOrb[i] == int(MO): # symmetry is hardcoded for now
            orbIndex = i
            break
      return orbIndex


# The path to the directory containing the folders
base_dir = '.'

with open('orbitals.csv', 'w', newline='') as csvfile:
    fieldnames = ['folder', 'distance', 'angle', 'overlap']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    # Loop over the folders in base_dir
    for folder_name in os.listdir(base_dir):
        # Continue only if the folder name ends with ".results" and does not contain "Region"
        if folder_name.endswith('.results') and "Region" not in folder_name:
            # Open the .rkf file
            kf = KFReader(os.path.join(base_dir, folder_name, 'adf.rkf'))

            # Extract the atom info and coordinates
            nAtoms, sAtoms, coordinates = GetAtom(kf)

            # Get coordinates of the atoms
            coord1 = coordinates[0]  # Atom 1
            coord11 = coordinates[10]  # Atom 11
            coord4 = coordinates[3]  # Atom 4

            # Calculate the distance and angle
            distance = calculate_distance(coord1, coord11)
            angle = calculate_angle(coord4, coord1, coord11)
            
            orb1_frag = 1
            orb1_MO = 17
            orb2_frag = 2
            orb2_MO = 7
            
            orb_index1 = GetOrbitalIndex(kf, orb1_frag, orb1_MO)  
            orb_index2 = GetOrbitalIndex(kf, orb2_frag, orb2_MO) 
            faIrrep = GetFaIrrep(kf)

            overlap = ReadOverlap(kf, faIrrep, orb_index1, orb_index2)

            # Write the data to the csv file
            writer.writerow({'folder': folder_name, 'distance': distance, 'angle': angle, 'overlap': overlap})


finish()
