import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# List of simulation IDs to compare
sim_ids = ["001", "002", "003", "004", "005"]  # Add more IDs as needed

# Load benchmark data
benchmark_data = pd.read_csv('MR_LBM/post/benchmark/ghia_ux_dy.csv', header=None)
y_d = benchmark_data.iloc[1:, 0].astype(float).values
u_x = benchmark_data.iloc[1:, 1:].astype(float).values.T
re_numbers = [f"Re= {int(value)}" for value in benchmark_data.iloc[0, 1:].values]

plt.figure(figsize=(10, 6))

# Plot benchmark data (e.g., for Re=10000 - adjust index as needed)
# benchmark index 0 - Re 100
# benchmark index 1 - Re 400
# benchmark index 2 - Re 1000
# benchmark index 3 - Re 3200
# benchmark index 4 - Re 5000
# benchmark index 5 - Re 7500
# benchmark index 6 - Re 10000

benchmark_index = 6  # Index for desired Reynolds number
plt.plot(u_x[benchmark_index], y_d, 'x', label=re_numbers[benchmark_index] + " (Benchmark)")

# Process and plot each simulation
for sim_id in sim_ids:
    # Paths to simulation files
    info_path = f"LDC/{sim_id}/info.txt"
    velocity_path = f"LDC/{sim_id}/velocity_x.bin"
    
    # Read parameters from info.txt
    nx_value = ny_value = umax_value = None
    try:
        with open(info_path, 'r') as f_info:
            for line in f_info:
                stripped = line.strip()
                
                if stripped.startswith('NX:'):
                    nx_value = int(stripped.split(':')[1].strip())
                elif stripped.startswith('NY:'):
                    ny_value = int(stripped.split(':')[1].strip())
                elif stripped.startswith('Umax:'):
                    umax_value = float(stripped.split(':')[1].strip())
    except FileNotFoundError:
        print(f"Info file not found for {sim_id}. Skipping.")
        continue
    
    # Validate parameters
    if None in (nx_value, ny_value, umax_value):
        print(f"Missing parameters in {info_path}. Skipping.")
        continue
    
    # Load and process velocity data
    try:
        with open(velocity_path, 'rb') as f_velocity:
            ux_sim = np.frombuffer(f_velocity.read(), dtype=np.float32)
    except FileNotFoundError:
        print(f"Velocity file not found for {sim_id}. Skipping.")
        continue
    
    # Generate normalized coordinates and plot
    y_sim = np.linspace(0, 1, ny_value)
    plt.plot(ux_sim / umax_value, y_sim, label=f"Simulation {sim_id}")

# Configure and show plot
plt.xlabel('u_x / Umax')
plt.ylabel('y / H')
plt.title('Velocity Profile Comparison')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()