import numpy as np
import matplotlib.pyplot as plt

# Define simulation IDs and line styles/colors
sim_ids = ['001', '002', '003']
colors = ['b', 'r', 'g']
linestyles = ['-', '--', '-.']

plt.figure(figsize=(10, 6))

for sim_id, color, ls in zip(sim_ids, colors, linestyles):
    # Construct file path for current simulation
    tke_path = f'LDC/{sim_id}/total_kinetic_energy.bin'
    
    try:
        # Read binary data
        with open(tke_path, 'rb') as f:
            data = np.frombuffer(f.read(), dtype=np.float32)
        
        # Extract time and kinetic energy
        t_star = data[::2]    # Even indices
        tke_sum = data[1::2]  # Odd indices
        
        # Plot with unique style and label
        plt.plot(t_star, tke_sum, 
                 color=color, 
                 linestyle=ls, 
                 linewidth=2,
                 label=f'Simulation {sim_id}')
    
    except FileNotFoundError:
        print(f"Warning: File not found - {tke_path}")

# Configure plot
plt.xlabel('Normalized Time (t*)')
plt.ylabel('Total Kinetic Energy')
plt.title('Kinetic Energy Evolution Comparison')
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()

# Save and show plot
##plt.savefig('kinetic_energy_comparison.png', dpi=300)
plt.show()