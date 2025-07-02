import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load and process data
print("Loading data...")

# Load data (assuming first column is y/d, subsequent columns are u_x at different Re/resolutions)
benchmark_data = pd.read_csv('MR_LBM/post/benchmark/ghia_ux_dy.csv', header=None)

velocity_path = "LDC/003/velocity_x.bin"

with open(velocity_path, 'rb') as f:
    sim_ux = np.frombuffer(f.read(), dtype=np.float32)

# Extract coordinates and velocities vou
y_d = benchmark_data.iloc[1:, 0].astype(float).values  # Normalized y-coordinate (y/Ny)
u_x = benchmark_data.iloc[1:, 1:].astype(float).values.T  # Transposed for easier plotting

# Generate y-coordinates from 0 to 1 with Ny points
y_sim = np.linspace(0, 1, 128)

# Get column headers (Re numbers or resolutions) if available
re_numbers = [f"Re= {int(value)}" for value in benchmark_data.iloc[0, 1:].values]

plt.plot(u_x[2], y_d, marker = 'x', linestyle= '', label = re_numbers[2])
plt.plot(sim_ux / 0.0256, y_sim)

plt.legend(title="Reynolds Number", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()