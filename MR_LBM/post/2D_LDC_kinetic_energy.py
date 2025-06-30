import numpy as np
import matplotlib.pyplot as plt

tke_path = 'LDC/001/total_kinetic_energy.bin'

# Read binary data
with open(tke_path, 'rb') as f:
    data = np.frombuffer(f.read(), dtype=np.float32)

# Reshape into pairs (t_star, tke_sum)
t_star = data[::2]    # Even indices: 0, 2, 4...
tke_sum = data[1::2]  # Odd indices: 1, 3, 5...

# Plot kinetic energy evolution
plt.figure(figsize=(10, 6))
plt.plot(t_star, tke_sum)
plt.xlabel('Normalized Time (t*)')
plt.ylabel('Total Kinetic Energy')
plt.title('Kinetic Energy Evolution')
plt.savefig('kinetic_energy.png', dpi=300)
plt.show()