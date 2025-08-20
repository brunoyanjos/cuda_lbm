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

coarse_rho_path = "GRID/001/coarse_rho.bin"
fine_rho_path = "GRID/001/fine_rho.bin"

coarse_ux_path = "GRID/001/coarse_ux.bin"
fine_ux_path = "GRID/001/fine_ux.bin"

coarse_uy_path = "GRID/001/coarse_uy.bin"
fine_uy_path = "GRID/001/fine_uy.bin"

coarse_mxx_path = "GRID/001/coarse_mxx.bin"
fine_mxx_path = "GRID/001/fine_mxx.bin"

coarse_mxy_path = "GRID/001/coarse_mxy.bin"
fine_mxy_path = "GRID/001/fine_mxy.bin"

coarse_myy_path = "GRID/001/coarse_myy.bin"
fine_myy_path = "GRID/001/fine_myy.bin"

with open(coarse_rho_path, 'rb') as f:
    coarse_rho = np.frombuffer(f.read(), dtype=np.float32)
    
with open(fine_rho_path, 'rb') as f:
    fine_rho = np.frombuffer(f.read(), dtype=np.float32)

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
    
calc_mxy_fine = (1 / ALPHA) * (coarse_mxy - coarse_ux * coarse_uy) + fine_ux[::2] * fine_uy[::2]
calc_mxy_coarse = ALPHA * (fine_mxy[::2] - fine_ux[::2] * fine_uy[::2]) + coarse_ux * coarse_uy
            
y_coarse = np.linspace(0, 1, NY_COARSE)
y_fine = np.linspace(0, 1, NY_FINE)

plt.plot(coarse_uy , y_coarse, label = "coarse")
plt.plot(fine_uy, y_fine, label = "fine")

plt.legend(title="grid kind", loc='lower right')

plt.tight_layout()

plt.show()
