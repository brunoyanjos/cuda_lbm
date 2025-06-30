import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load and process data
print("Loading data...")

# Load data (assuming first column is y/d, subsequent columns are u_x at different Re/resolutions)
data = pd.read_csv('MR_LBM/post/benchmark/ghia_ux_dy.csv', header=None)

# Extract coordinates and velocities vou
y_d = data.iloc[1:, 0].astype(float).values  # Normalized y-coordinate (y/Ny)
u_x = data.iloc[1:, 1:].astype(float).values.T  # Transposed for easier plotting

# Get column headers (Re numbers or resolutions) if available
re_numbers = [f"Re= {int(value)}" for value in data.iloc[0, 1:].values]

plt.plot(u_x[0], y_d, marker = 'x', linestyle= '', label = re_numbers[0])
plt.plot(u_x[1], y_d, marker = '.', linestyle= '', label = re_numbers[1])

plt.legend(title="Reynolds Number", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()