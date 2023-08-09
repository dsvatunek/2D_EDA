"""
ADF Input File Generator
Author: Dennis Svatunek

Description:
Script to produce ADF 2023 input files for EDA analysis from multi-structure XYZ file: structures.xyz
"""

import os


def read_multi_xyz(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    molecules = []
    i = 0
    while i < len(lines):
        num_atoms = int(lines[i].strip())
        molecule_lines = lines[i + 2:i + 2 + num_atoms]
        atoms = [line.split() for line in molecule_lines]
        molecules.append(atoms)
        i += num_atoms + 2

    return molecules


def write_adf_input(molecule_idx, atoms):

    # Writing the main ADF input file
    with open(f"{molecule_idx}.run", 'w') as file:
        # Writing the dependencies
        file.write("#!/bin/sh\n\n")
        file.write(f"# dependency: {molecule_idx}.Region_1 {molecule_idx}.Region_1.results/adf.rkf Region_1.rkf\n")
        file.write(f"# dependency: {molecule_idx}.Region_2 {molecule_idx}.Region_2.results/adf.rkf Region_2.rkf\n\n")
        
        file.write(f'"$AMSBIN/ams" << eor\n\n')
        file.write(f"Task SinglePoint\n")
        file.write(f"System\n")
        file.write(f"    Atoms\n")
        
        # Atoms for fragment 1
        for i in range(10):
            file.write(f"        {atoms[i][0]} {' '.join(atoms[i][1:])} region=Region_1 adf.f=Region_1\n")
        # Atoms for fragment 2
        for i in range(10, 12):
            file.write(f"        {atoms[i][0]} {' '.join(atoms[i][1:])} region=Region_2 adf.f=Region_2\n")
        
        file.write(f"End\n")
        file.write(f"    Charge -1\n")
        file.write(f"End\n")
        file.write(f"\n")
        file.write(f"Engine ADF\n")
        file.write(f"    Basis\n")
        file.write(f"        Type TZ2P\n")
        file.write(f"        Core None\n")
        file.write(f"    End\n")
        file.write(f"    Fragments\n")
        file.write(f"      Region_1 Region_1.rkf \n")
        file.write(f"      Region_2 Region_2.rkf \n")
        file.write(f"    End\n")
        file.write(f"    XC\n")
        file.write(f"        MetaHybrid M06-2X\n")
        file.write(f"    End\n")
        file.write(f"    Symmetry NOSYM\n")
        file.write(f"    NumericalQuality VeryGood\n")
        file.write(f"EndEngine\n")
        file.write(f"eor\n")
    
    # Writing the fragment 1 input file
    with open(f"{molecule_idx}.region_1.run", 'w') as file:
        file.write("#!/bin/sh\n\n")
        file.write(f'"$AMSBIN/ams" << eor\n\n')
        file.write(f"Task SinglePoint\n")
        file.write(f"System\n")
        file.write(f"    Atoms\n")
        
        # Atoms for fragment 1
        for i in range(10):
            file.write(f"        {atoms[i][0]} {' '.join(atoms[i][1:])} region=Region_1\n")
        
        file.write(f"End\n")
        file.write(f"    Symmetrize Yes\n")
        file.write(f"End\n")
        file.write(f"Symmetry\n")
        file.write(f"    SymmetrizeTolerance 1.0e-7\n")
        file.write(f"End\n")
        file.write(f"\n")
        file.write(f"Engine ADF\n")
        file.write(f"    Basis\n")
        file.write(f"        Type TZ2P\n")
        file.write(f"        Core None\n")
        file.write(f"    End\n")
        file.write(f"    XC\n")
        file.write(f"        MetaHybrid M06-2X\n")
        file.write(f"    End\n")
        file.write(f"    Title Region_1 fragment\n")
        file.write(f"    Symmetry NOSYM\n")
        file.write(f"    NumericalQuality VeryGood\n")
        file.write(f"EndEngine\n")
        file.write(f"eor\n")
    
    # Writing the fragment 2 input file
    with open(f"{molecule_idx}.region_2.run", 'w') as file:
        file.write("#!/bin/sh\n\n")
        file.write(f'"$AMSBIN/ams" << eor\n\n')
        file.write(f"Task SinglePoint\n")
        file.write(f"System\n")
        file.write(f"    Atoms\n")
        
        # Atoms for fragment 2
        for i in range(10, 12):
            file.write(f"        {atoms[i][0]} {' '.join(atoms[i][1:])} region=Region_2\n")
        
        file.write(f"End\n")
        file.write(f"    Symmetrize Yes\n")
        file.write(f"    Charge -1\n")
        file.write(f"End\n")
        file.write(f"Symmetry\n")
        file.write(f"    SymmetrizeTolerance 1.0e-7\n")
        file.write(f"End\n")
        file.write(f"\n")
        file.write(f"Engine ADF\n")
        file.write(f"    Basis\n")
        file.write(f"        Type TZ2P\n")
        file.write(f"        Core None\n")
        file.write(f"    End\n")
        file.write(f"    XC\n")
        file.write(f"        MetaHybrid M06-2X\n")
        file.write(f"    End\n")
        file.write(f"    Title Region_2 fragment\n")
        file.write(f"    Symmetry NOSYM\n")
        file.write(f"    NumericalQuality VeryGood\n")
        file.write(f"EndEngine\n")
        file.write(f"eor\n")


# Reading the multi xyz file
molecules = read_multi_xyz("structures.xyz")

# Creating ADF input files for each molecule in the multi xyz file
for idx, molecule in enumerate(molecules):
    write_adf_input(idx, molecule)
