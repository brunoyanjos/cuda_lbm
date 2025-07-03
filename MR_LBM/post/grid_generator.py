import numpy as np
import matplotlib.pyplot as plt

Nx = 5
Ny = Nx
Ni = 1
Ne = 1

# Define grid parameters
coarse_step = 1.0  # Coarse grid spacing
fine_step = 1/2    # Fine grid spacing

# Generate coarse grid points

coarse_x, coarse_y = np.meshgrid(
    np.arange(Ni + Ne, (Nx + Ni + Ne), coarse_step),
    np.arange(Ni + Ne, (Ny + Ni + Ne), coarse_step)
)

# Generate fine grid points
fine_x, fine_y = np.meshgrid(
    np.arange(0, (Nx + 2 * Ni + 2 * Ne) - fine_step, fine_step),
    np.arange(0, (Ny + 2 * Ni + 2 * Ne) - fine_step, fine_step)
)

# Create plot
plt.figure(figsize=(10, 8))

# Plot fine grid (smaller dots)
plt.scatter(fine_x, fine_y, s=20, color='red', label='Fine Grid')

# Plot coarse grid (larger dots)
plt.scatter(coarse_x, coarse_y, s=80, color='black', label='Coarse Grid')

# Configure plot appearance
plt.title('Overlaid Coarse and Fine Grids', fontsize=14)
plt.xlabel('X-axis', fontsize=12)
plt.ylabel('Y-axis', fontsize=12)
plt.xticks(np.arange(0, (Nx + 2 * Ni + 2 * Ne), 1))
plt.yticks(np.arange(0, (Ny + 2 * Ni + 2 * Ne), 1))
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend(loc='upper right')
plt.gca().set_aspect('equal')  # Ensure equal aspect ratio

# Show the plot
plt.tight_layout()
plt.show()