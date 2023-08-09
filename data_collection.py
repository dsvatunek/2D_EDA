"""
EDA data collector
Author: Dennis Svatunek
Date: 2023-08-09

Description:
This script reads EDA energy components from the EDA results and compiles into CSV.
"""

import glob
import re
import math

fragment1_energy = -1781.02  # Acetone energy isolated relaxed structure
fragment2_energy = -562.35   # Cyanide energy isolated relaxed structure

energy_unit = 2  # 2 for kcal/mol


def find_value_in_file(file_path, target_string, index):
    with open(file_path, 'r') as file:
        for line in file:
            if target_string in line:
                # Extract all numbers from the line
                numbers = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                # Return the number at the specified index
                return float(numbers[index])


def extract_geometry(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Find start and end of coordinates block
    start = None
    end = None
    for i, line in enumerate(lines):
        if 'Atoms' in line:
            start = i + 1
        if 'End' in line:
            end = i
            break

    # Extract coordinates
    coords = [line.split() for line in lines[start:end]]
    atom1 = [float(x) for x in coords[0][1:4]]
    atom11 = [float(x) for x in coords[10][1:4]]
    atom2 = [float(x) for x in coords[1][1:4]]

    # Calculate distance 1/11
    distance = math.sqrt((atom1[0] - atom11[0])**2 + (atom1[1] - atom11[1])**2 + (atom1[2] - atom11[2])**2)

    # Calculate angle 2/1/11
    a = [atom2[i] - atom1[i] for i in range(3)]
    b = [atom11[i] - atom1[i] for i in range(3)]
    dot_product = sum(a[i] * b[i] for i in range(3))
    mag_a = math.sqrt(sum(x**2 for x in a))
    mag_b = math.sqrt(sum(x**2 for x in b))
    angle = math.acos(dot_product / (mag_a * mag_b))

    return distance, math.degrees(angle)


with open('results.csv', 'w') as results:
    results.write('Index,Distance 1-11,Angle 2-1-11,Total Pauli Repulsion,Electrostatic Interaction,Total Orbital Interactions,Total Interaction Energy,Strain Fragment1,Strain Fragment2,Total Strain,Relative Energy\n')

    for main_file in glob.glob('[0-9]*.out'):
        if '.Region_1.' in main_file or '.Region_2.' in main_file:
            continue

        index = main_file.split('.')[0]

        run_file = main_file.replace('.out', '.run')
        distance, angle = extract_geometry(run_file)

        total_pauli_repulsion = find_value_in_file(main_file, "Total Pauli Repulsion:", energy_unit)
        electrostatic_interaction = find_value_in_file(main_file, "Electrostatic Interaction:", energy_unit)
        total_orbital_interactions = find_value_in_file(main_file, "Total Orbital Interactions:", energy_unit)
        total_interaction_energy = find_value_in_file(main_file, "Total Bonding Energy:", energy_unit)

        # Read fragment files
        fragment1_file = main_file.replace('.out', '.Region_1.out')
        fragment2_file = main_file.replace('.out', '.Region_2.out')

        strain_fragment1 = find_value_in_file(fragment1_file, "Total Bonding Energy:", energy_unit) - fragment1_energy
        strain_fragment2 = find_value_in_file(fragment2_file, "Total Bonding Energy:", energy_unit) - fragment2_energy
        total_strain = strain_fragment1 + strain_fragment2
        sum_interaction_and_total_strain = total_interaction_energy + total_strain

        results.write(f"{index},{distance},{angle},{total_pauli_repulsion},{electrostatic_interaction},{total_orbital_interactions},{total_interaction_energy},{strain_fragment1},{strain_fragment2},{total_strain},{sum_interaction_and_total_strain}\n")
