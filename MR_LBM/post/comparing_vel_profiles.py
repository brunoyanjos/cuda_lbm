import numpy as np
import matplotlib.pyplot as plt

RE = 1000

GRID_RATIO = 2

U_MAX = 0.0256

NY_COARSE = 33
NY_FINE = NY_COARSE * GRID_RATIO

VISC_FINE = U_MAX * (NY_FINE - 1) / RE
VISC_COARSE = U_MAX * (NY_COARSE - 1) / RE

TAU_FINE = 0.5 + 3.0 * VISC_FINE
TAU_COARSE = 0.5 + 3.0 * VISC_COARSE

OMEGA_FINE = 1.0 / TAU_FINE
OMEGA_COARSE = 1.0 / TAU_COARSE

ALPHA = OMEGA_COARSE / OMEGA_FINE * GRID_RATIO

coarse_ux_path = "GRID/001/coarse_ux.bin"
fine_ux_path = "GRID/001/fine_ux.bin"

coarse_uy_path = "GRID/001/coarse_uy.bin"
fine_uy_path = "GRID/001/fine_uy.bin"

coarse_mxy_path = "GRID/001/coarse_mxy.bin"
fine_mxy_path = "GRID/001/fine_mxy.bin"

with open(coarse_ux_path, 'rb') as f:
    coarse_ux = np.frombuffer(f.read(), dtype=np.float32)
    
with open(fine_ux_path, 'rb') as f:
    fine_ux = np.frombuffer(f.read(), dtype=np.float32)
    
with open(coarse_uy_path, 'rb') as f:
    coarse_uy = np.frombuffer(f.read(), dtype=np.float32)
    
with open(fine_uy_path, 'rb') as f:
    fine_uy = np.frombuffer(f.read(), dtype=np.float32)
    
with open(coarse_mxy_path, 'rb') as f:
    coarse_mxy = np.frombuffer(f.read(), dtype=np.float32)
    
with open(fine_mxy_path, 'rb') as f:
    fine_mxy = np.frombuffer(f.read(), dtype=np.float32)
    
calc_mxy = ALPHA * (coarse_mxy - coarse_ux * coarse_uy) + fine_ux[::2] * fine_uy[::2]
    
print(fine_ux)
print(fine_ux[::2])
            
y_coarse = np.linspace(0, 1, NY_COARSE)
y_fine = np.linspace(0, 1, NY_FINE)

plt.plot(coarse_mxy , y_coarse, label = "coarse")
plt.plot(fine_mxy , y_fine, label = "fine")
plt.plot(calc_mxy , y_coarse, label = "calc_fine")

plt.legend(title="grid kind", loc='lower right')

plt.tight_layout()

plt.show()
