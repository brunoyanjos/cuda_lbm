import numpy as np
import matplotlib.pyplot as plt

# === Estilo LaTeX (opcional) ===
plt.rcParams.update({
    "text.usetex": False,  # usa o motor interno, não o LaTeX real
    "font.family": "serif",
    "mathtext.fontset": "cm",  # usa fontes Computer Modern
})

RE = 100

GRID_RATIO = 2

U_MAX = 0.01

NY_COARSE = 64
NY_FINE = NY_COARSE * GRID_RATIO - 1

VISC_FINE = U_MAX * (NY_FINE - 1) / RE
VISC_COARSE = U_MAX * (NY_COARSE - 1) / RE

TAU_FINE = 0.5 + 3.0 * VISC_FINE
TAU_COARSE = 0.5 + 3.0 * VISC_COARSE

OMEGA_FINE = 1.0 / TAU_FINE
OMEGA_COARSE = 1.0 / TAU_COARSE

ALPHA = OMEGA_FINE * GRID_RATIO / OMEGA_COARSE

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
    
with open(coarse_mxx_path, 'rb') as f:
    coarse_mxx = np.frombuffer(f.read(), dtype=np.float32)
    
with open(fine_mxx_path, 'rb') as f:
    fine_mxx = np.frombuffer(f.read(), dtype=np.float32)    

with open(coarse_mxy_path, 'rb') as f:
    coarse_mxy = np.frombuffer(f.read(), dtype=np.float32)
    
with open(fine_mxy_path, 'rb') as f:
    fine_mxy = np.frombuffer(f.read(), dtype=np.float32)
    
with open(coarse_myy_path, 'rb') as f:
    coarse_myy = np.frombuffer(f.read(), dtype=np.float32)
    
with open(fine_myy_path, 'rb') as f:
    fine_myy = np.frombuffer(f.read(), dtype=np.float32)
  
y_coarse = np.linspace(0, 1, NY_COARSE)
y_fine = np.linspace(0, 1, NY_FINE)


green = '#149B55'
dark_green = '#054B28'

red = "#F04137"
coral = '#FF6968'

# -------------------------------------------------------------------------
# ------------------------------ RHO PROFILE ------------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

plt.plot(coarse_rho, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(fine_rho, y_fine, label="Fine grid", 
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$\rho$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_rho_normalizado.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------------ UX PROFILE -------------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

u_norm = np.max([np.max(coarse_ux), np.max(fine_ux)]) 

plt.plot(coarse_ux / u_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(fine_ux / u_norm, y_fine, label="Fine grid", 
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$u_x / U_x^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_ux_normalizado.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------------ UY PROFILE -------------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

u_norm = np.max([np.max(abs(coarse_uy)), np.max(abs(fine_uy))]) 

plt.plot(coarse_uy / u_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(fine_uy / u_norm, y_fine, label="Fine grid", 
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$u_y / U_y^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_uy_normalizado.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------------ MXX PROFILE ------------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

m_norm = np.max([np.max(abs(coarse_mxx)), np.max(abs(fine_mxx))]) 

plt.plot(coarse_mxx / m_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(fine_mxx / m_norm, y_fine, label="Fine grid", 
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$m_{xx} / m_{xx}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_mxx_normalizado.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------------ MXY PROFILE ------------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

m_norm = np.max([np.max(abs(coarse_mxy)), np.max(abs(fine_mxy))]) 

plt.plot(coarse_mxy / m_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(fine_mxy / m_norm, y_fine, label="Fine grid", 
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$m_{xy} / m_{xy}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_mxy_normalizado.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------------ MYY PROFILE ------------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

m_norm = np.max([np.max(abs(coarse_myy)), np.max(abs(fine_myy))]) 

plt.plot(coarse_myy / m_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(fine_myy / m_norm, y_fine, label="Fine grid", 
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$m_{yy} / m_{yy}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_myy_normalizado.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------ MXX PROFILE CORRECTION -------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

calc_mxx_coarse = ALPHA * (fine_mxx[::2] - fine_ux[::2] * fine_ux[::2]) + coarse_ux * coarse_ux

m_norm = np.max([np.max(abs(coarse_mxx)), np.max(abs(fine_mxx))]) 

plt.plot(coarse_mxx / m_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(calc_mxx_coarse / m_norm, y_coarse, label="Evaluate Coarse", 
        linestyle='--', linewidth=1.8, color=red)

plt.xlabel(r"$m_{xx} / m_{xx}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_mxx_coarse_evaluate.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

calc_mxx_fine = (1 / ALPHA) * (coarse_mxx - coarse_ux * coarse_ux) + fine_ux[::2] * fine_ux[::2]

m_norm = np.max([np.max(abs(calc_mxx_fine)), np.max(abs(fine_mxx))])

plt.plot(calc_mxx_fine / m_norm, y_coarse, label="Evaluate Fine", 
         marker='x', markersize=5, linewidth=1.6, color=dark_green)
plt.plot(fine_mxx / m_norm, y_fine, label="Fine grid",
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$m_{xx} / m_{xx}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_mxx_fine_evaluate.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------ MXY PROFILE CORRECTION -------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

calc_mxy_coarse = ALPHA * (fine_mxy[::2] - fine_ux[::2] * fine_uy[::2]) + coarse_ux * coarse_uy

m_norm = np.max([np.max(abs(coarse_mxy)), np.max(abs(calc_mxy_coarse))])

plt.plot(coarse_mxy / m_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(calc_mxy_coarse / m_norm, y_coarse, label="Evaluate Coarse", 
        linestyle='--', linewidth=1.8, color=red)

plt.xlabel(r"$m_{xy} / m_{xy}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_mxy_coarse_evaluate.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

calc_mxy_fine = (1 / ALPHA) * (coarse_mxy - coarse_ux * coarse_uy) + fine_ux[::2] * fine_uy[::2]

m_norm = np.max([np.max(abs(calc_mxy_fine)), np.max(abs(fine_mxy))])

plt.plot(calc_mxy_fine / m_norm, y_coarse, label="Evaluate Fine", 
         marker='x', markersize=5, linewidth=1.6, color=dark_green)
plt.plot(fine_mxy / m_norm, y_fine, label="Fine grid",
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$m_{xy} / m_{xy}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_mxy_fine_evaluate.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------
# ------------------------ MYY PROFILE CORRECTION -------------------------
# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

calc_myy_coarse = ALPHA * (fine_myy[::2] - fine_uy[::2] * fine_uy[::2]) + coarse_uy * coarse_uy

m_norm = np.max([np.max(abs(coarse_myy)), np.max(abs(calc_myy_coarse))])

plt.plot(coarse_myy / m_norm, y_coarse, label="Coarse grid", 
         marker='x', markersize=5, linewidth=1.6, color=green)
plt.plot(calc_myy_coarse / m_norm, y_coarse, label="Evaluate Coarse", 
        linestyle='--', linewidth=1.8, color=red)

plt.xlabel(r"$m_{yy} / m_{yy}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_myy_coarse_evaluate.png", dpi=300, bbox_inches='tight')
plt.close()

# -------------------------------------------------------------------------

plt.figure(figsize=(6, 4), dpi=300)

calc_myy_fine = (1 / ALPHA) * (coarse_myy - coarse_uy * coarse_uy) + fine_uy[::2] * fine_uy[::2]

m_norm = np.max([np.max(abs(calc_myy_fine)), np.max(abs(fine_myy))])

plt.plot(calc_myy_fine / m_norm, y_coarse, label="Evaluate Fine", 
         marker='x', markersize=5, linewidth=1.6, color=dark_green)
plt.plot(fine_myy/ m_norm, y_fine, label="Fine grid",
         linestyle='--', linewidth=1.8, color=coral)

plt.xlabel(r"$m_{yy} / m_{yy}^{\mathrm{max}}$", fontsize=12)
plt.ylabel(r"$y / L$", fontsize=12)

plt.legend(fontsize=10, title_fontsize=10, loc='lower right')
plt.grid(True, linestyle=':', linewidth=0.8, alpha=0.7)

plt.tight_layout()
plt.savefig("perfil_myy_fine_evaluate.png", dpi=300, bbox_inches='tight')
plt.close()
