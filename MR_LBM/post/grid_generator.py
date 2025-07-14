import numpy as np
import matplotlib.pyplot as plt

Nx = 32
Ny = Nx
Ni = 1
Ne = 3

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

Nx_fine_grid = (Nx_total_size - 1) * fine_to_coarse + 1
Ny_fine_grid = (Ny_total_size - 1) * fine_to_coarse + 1

Nx_coarse_grid = Nx + 2 * Ni
Ny_coarse_grid = Ny + 2 * Ni

Nx_fine_width = (Ni + Ne) * fine_to_coarse + 1
Ny_fine_width = (Ni + Ne) * fine_to_coarse + 1

print(f'Nx_fine_grid: {Nx_fine_grid}, Ny_fine_grid: {Ny_fine_grid}')
print(f'Nx_coarse_grid: {Nx_coarse_grid}, Ny_coarse_grid: {Ny_coarse_grid}')
print(f'Nx_fine_width: {Nx_fine_width}, Ny_fine_width: {Ny_fine_width}')

for i in range(Nx_fine_grid):    
    for j in range(Ny_fine_grid):
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
            
total_coarse_points = len(coarse_grid_x) + len(intersection_grid_x)
total_fine_points = len(fine_grid_x)
total_points = total_coarse_points + total_fine_points
similarity = (total_points) ** (1/2)
current_accuracy = Nx_fine_grid * Ny_fine_grid

print(f'total_coarse_points: {total_coarse_points}')
print(f'total_fine_points: {2 * (Ny_fine_width * Nx_fine_grid + Nx_fine_width * (Ny_fine_grid - 2 * Ny_fine_width))}')
print(f'total_fine_points: {total_fine_points}')
print(f'similarity: {similarity}')
print(f'savings: {((current_accuracy - total_points)/current_accuracy * 100):.2f}%')
        
          
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
plt.xticks(np.arange(0, (Nx + 2 * Ni + 2 * Ne)))
plt.yticks(np.arange(0, (Ny + 2 * Ni + 2 * Ne)))
plt.grid(True, linestyle='--', alpha=0.3)
plt.gca().set_aspect('equal')  # Ensure equal aspect ratio

plt.tight_layout()
plt.show()