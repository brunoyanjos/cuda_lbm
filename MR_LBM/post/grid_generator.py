import numpy as np
import matplotlib.pyplot as plt

Nx = 3
Ny = Nx
Ni = 1
Ne = 1

# Define grid parameters
coarse_step = 1.0  # Coarse grid spacing
fine_to_coarse = 2
fine_step = 1 / fine_to_coarse  # Fine grid spacing

coarse_grid_x = []
coarse_grid_y = []

intersection_grid_x = []
intersection_grid_y = []

fine_grid_x = []
fine_grid_y = []

Nx_total_size = (Nx + 2 * Ni + 2 * Ne)
Ny_total_size = (Ny + 2 * Ni + 2 * Ne)

start_of_coarse = Ni + Ne
stop_coarse = start_of_coarse + Nx - 1

print(f"start_of_coarse: {start_of_coarse}, stop_coarse: {stop_coarse}")

for i in range((Nx_total_size - 1) * fine_to_coarse + 1):
    for j in range((Ny_total_size - 1) * fine_to_coarse + 1):
        x = i * fine_step
        y = j * fine_step
        
        if x % 1 == 0 and y % 1 == 0: 
            if x >= start_of_coarse and x <= stop_coarse and y >= start_of_coarse and y <= stop_coarse:
                coarse_grid_x.append(x)
                coarse_grid_y.append(y)
                
            elif x >= start_of_coarse - Ni and x <= stop_coarse + Ni and y >= start_of_coarse - Ni and y <= stop_coarse + Ni:
                intersection_grid_x.append(x)
                intersection_grid_y.append(y)
             
        if x <= start_of_coarse or x >= stop_coarse or y <= start_of_coarse or y >= stop_coarse:
            fine_grid_x.append(x)
            fine_grid_y.append(y)
        
          
plt.scatter(fine_grid_x, fine_grid_y, color='red', s= 10)
plt.scatter(coarse_grid_x, coarse_grid_y, color='black', s= 40)
plt.scatter(intersection_grid_x, intersection_grid_y, s=40,                  # Marker size
            facecolors='none',       # Transparent center
            edgecolors='black',        # Border color
            linewidths=1,           # Border thickness
            marker='o')              # Circle shape
plt.xlabel('X-axis', fontsize=12)
plt.ylabel('Y-axis', fontsize=12)
plt.xticks(np.arange(0, (Nx + 2 * Ni + 2 * Ne)))
plt.yticks(np.arange(0, (Ny + 2 * Ni + 2 * Ne)))
plt.grid(True, linestyle='--', alpha=0.3)
plt.gca().set_aspect('equal')  # Ensure equal aspect ratio

plt.tight_layout()
plt.show()