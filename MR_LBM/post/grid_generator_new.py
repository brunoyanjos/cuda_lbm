import numpy as np
import matplotlib.pyplot as plt

Nx = 9
Ny = Nx
Ni = 1

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

Nx_coarse = int(Nx / 2) + 1
Ny_coarse = Ny

Nx_fine = Nx
Ny_fine = Ny_coarse * fine_to_coarse - 1

for i in range(Nx * 2 - 1):
   for j in range(Ny * 2 - 1):
       x = i * 0.5
       y = j * 0.5
       
       if x % 1 == 0 and y % 1 == 0 and x < Nx_coarse:
        coarse_grid_x.append(x)
        coarse_grid_y.append(y)
       elif x + 1 >= Nx_coarse:
        fine_grid_x.append(x)
        fine_grid_y.append(y)
               
          
plt.scatter(fine_grid_x, fine_grid_y, color='red', s= 10)
plt.scatter(coarse_grid_x, coarse_grid_y, color='black', s= 40)
plt.scatter(intersection_grid_x, intersection_grid_y, 
            s=40,                    # Marker size
            facecolors='none',       # Transparent center
            edgecolors='black',      # Border color
            linewidths=1,            # Border thickness
            marker='o')              # Circle shape
plt.xlabel('X-axis', fontsize=12)
plt.ylabel('Y-axis', fontsize=12)
plt.xticks(np.arange(0, Nx))
plt.yticks(np.arange(0, Ny))
plt.grid(True, linestyle='--', alpha=0.3)
plt.gca().set_aspect('equal')  # Ensure equal aspect ratio

plt.tight_layout()
plt.show()