import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Data import (change if needed)
fine_data = np.loadtxt('grid_fine001.dat')
fine_data = np.delete(fine_data, 0, axis=1)
fine_dataT = fine_data.T

coarse_data = np.loadtxt('grid_coarse001.dat')
coarse_data = np.delete(coarse_data, 0, axis=1)
coarse_dataT = coarse_data.T

Nx = 5
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

# Graphs pré-config
vmin = min(coarse_data.min(), fine_data.min())
vmax = max(coarse_data.max(), fine_data.max())

# Normalizing
xmax_real = (Nx * 2 - 2) * 0.5  # Biggest X value
ymax_real = (Ny * 2 - 2) * 0.5

#%%
for i in range(Nx * 2 - 1):
   for j in range(Ny * 2 - 1):
       
       x = i * 0.5
       y = j * 0.5
       
       x_norm = x / xmax_real
       y_norm = y / ymax_real
       
       if x % 1 == 0 and y % 1 == 0 and x < Nx_coarse:
           
        coarse_grid_x.append(x_norm)
        coarse_grid_y.append(y_norm)

        
       if x + 1 >= Nx_coarse:
           
        fine_grid_x.append(x_norm)
        fine_grid_y.append(y_norm)
        

#%% DOT GRAPHIC
           
fig, ax1 = plt.subplots(figsize = (6,6), dpi=300)          

plt.scatter(fine_grid_x, fine_grid_y, c=fine_dataT, cmap='viridis', s=50, vmin=vmin, vmax=vmax)
plt.scatter(coarse_grid_x, coarse_grid_y, c=coarse_dataT, cmap='viridis', s=110, vmin=vmin, vmax=vmax)
plt.colorbar(shrink=0.9)

plt.tight_layout()
plt.show()

#%% CONTOURF GRAPHIC

coarse_x_unique = np.unique(coarse_grid_x)
coarse_y_unique = np.unique(coarse_grid_y)
Xc, Yc = np.meshgrid(coarse_x_unique, coarse_y_unique)

    
fine_x_unique = np.unique(fine_grid_x)
fine_y_unique = np.unique(fine_grid_y)
Xf, Yf = np.meshgrid(fine_x_unique, fine_y_unique)

    
fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

# Contourf
cf1 = ax.contourf(Xc, Yc, coarse_data, levels=50, cmap='viridis', vmin=vmin, vmax=vmax, corner_mask=True)
cf2 = ax.contourf(Xf, Yf, fine_data, levels=50, cmap='viridis', vmin=vmin, vmax=vmax, corner_mask=True)

plt.colorbar(cf2, ax=ax, shrink=0.7)

# Dots
plt.scatter(fine_grid_x, fine_grid_y, color='red', s= 10, alpha=0.6)
plt.scatter(coarse_grid_x, coarse_grid_y, color='black', s= 40,alpha=0.6)


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_aspect('equal')
# ax.set_title('Malhas sobrepostas com colormap')
ax.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.show()

#%% DOTS GRAPHIC (NO VALUES)

fig, ax1 = plt.subplots(figsize = (4,4), dpi=300)  

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
# plt.xticks(np.arange(0, Nx))
# plt.yticks(np.arange(0, Ny))
plt.grid(True, linestyle='--', alpha=0.3)
plt.gca().set_aspect('equal')  # Ensure equal aspect ratio

plt.tight_layout()
plt.show()
