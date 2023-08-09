"""
ADF print cartesian
Author: Dennis Svatunek

Description:
This script prints EDA energy components on a 2D surface in cartesian coordinates
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy.interpolate import griddata
import os


def create_plot(x, y, z_values, title, filename):
    # Create grid and interpolate data
    grid_x, grid_y = np.mgrid[min(min(x), 0):max(x):100j, min(min(y), 0):max(y):100j]
    grid_z = griddata((x, y), z_values, (grid_x, grid_y), method='cubic')

    # Plot
    plt.figure(figsize=(8, 6))
    num_levels = 500  # Increase this for smoother transitions
    contourf = plt.contourf(grid_x, grid_y, grid_z, levels=num_levels, cmap='rainbow')
    plt.colorbar(contourf, label='Energy')
    plt.contour(grid_x, grid_y, grid_z, levels=40, colors='grey', linewidths=0.5)  # Add contour lines
    plt.xlabel('x / Angstrom')
    plt.ylabel('y / Angstrom')
    plt.title(title)

    # Plot origin
    plt.scatter(0, 0, color='red', marker='o', label='Origin (0,0)')

    # Add lines originating from (0,0) at angles of 155, 111, and 80 degrees
    line_length = max(max(x), max(y))
    x_155 = line_length * np.cos(np.radians(155))
    y_155 = line_length * np.sin(np.radians(155))
    x_111 = line_length * np.cos(np.radians(111.6))
    y_111 = line_length * np.sin(np.radians(111.6))
    x_70 = line_length * np.cos(np.radians(70))
    y_70 = line_length * np.sin(np.radians(70))

    plt.plot([0, x_155], [0, y_155], 'k-', label='155 degrees', linewidth=0.5)
    plt.plot([0, x_111], [0, y_111], 'b--', label='111 degrees')
    plt.plot([0, x_70], [0, y_70], 'k-', label='70 degrees', linewidth=0.5)
    
    # for radius in [1.5, 1.7, 1.9]:
        # arc = patches.Arc((0, 0), 2*radius, 2*radius, angle=0, theta1=70, theta2=170, color='black', linestyle='dashed')
        # plt.gca().add_patch(arc)
    # plt.legend()
    
    # Set the x and y scales to be equal
    plt.axis('equal')

    # Save plot to PNG file with higher resolution
    plt.savefig(os.path.join('graphs', filename), dpi=300)
    plt.close()


# Create directory for graphs if it doesn't exist
if not os.path.exists('graphs'):
    os.makedirs('graphs')

# Read data from the CSV file
data = []
with open('results.csv', 'r') as file:
    header = next(file).strip().split(',')
    for line in file:
        data.append(list(map(float, line.strip().split(','))))

# Unpack data
data = np.array(data)
distances = data[:, 1]
angles = data[:, 2]

# Convert polar to Cartesian coordinates for better plotting
angles_radians = np.radians(angles)
x = distances * np.cos(angles_radians)
y = distances * np.sin(angles_radians)

# Loop through energy components and create plots
for i in range(3, len(header)):
    energy_values = data[:, i]
    title = header[i].replace('_', ' ').title()
    create_plot(x, y, energy_values, title, f'{title.replace(" ", "_")}.png')
