import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ==================================================
# CONFIGURAÇÕES GERAIS
# ==================================================
sim_id = "004"
base_path = Path("ANNUL") / sim_id

path_uy_0 = base_path / "uy_theta_0.bin"
path_ux_90 = base_path / "ux_theta_90.bin"
path_uy_180 = base_path / "uy_theta_180.bin"
path_ux_270 = base_path / "ux_theta_270.bin"

out_dir = Path("output") / sim_id
out_dir.mkdir(parents=True, exist_ok=True)

dtype = np.float32

# ==================================================
# PARÂMETROS FÍSICOS
# ==================================================
Ri = 64.0
Ro = 128.0

Omega_i = 0.0256 / Ri
Omega_o = 0.0

# ==================================================
# LEITURA DOS DADOS (UM PERFIL POR ARQUIVO)
# ==================================================
u_0 = np.fromfile(path_uy_0, dtype=dtype)
u_90 = np.fromfile(path_ux_90, dtype=dtype)
u_180 = np.fromfile(path_uy_180, dtype=dtype)
u_270 = np.fromfile(path_ux_270, dtype=dtype)

# Checagem básica
assert (
    len(u_0) == len(u_90) == len(u_180) == len(u_270)
), "Perfis com tamanhos diferentes!"

# ==================================================
# CORREÇÃO DE ORIENTAÇÃO RADIAL
# ==================================================
u_0 = u_0[::-1]
u_270 = u_270[::-1]

# ==================================================
# RAIO
# ==================================================
N = len(u_0)
r = np.linspace(Ri, Ro, N)

# ==================================================
# SOLUÇÃO ANALÍTICA
# ==================================================
A = (Omega_o * Ro**2 - Omega_i * Ri**2) / (Ro**2 - Ri**2)
B = ((Omega_i - Omega_o) * Ri**2 * Ro**2) / (Ro**2 - Ri**2)
u_ana = A * r + B / r


# ==================================================
# ERROS
# ==================================================
def rms(x):
    return np.sqrt(np.mean(x**2))


err_0 = u_0 - u_ana
err_90 = u_90 - u_ana
err_180 = u_180 - u_ana
err_270 = u_270 - u_ana

# ==================================================
# PLOT — LBM × ANALÍTICO
# ==================================================
fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)

profiles = [
    (axs[0, 0], r"$\theta = 0^\circ$", u_0, err_0),
    (axs[0, 1], r"$\theta = 90^\circ$", u_90, err_90),
    (axs[1, 0], r"$\theta = 180^\circ$", u_180, err_180),
    (axs[1, 1], r"$\theta = 270^\circ$", u_270, err_270),
]

for ax, title, u_num, err in profiles:
    ax.plot(r, u_num, lw=2, label="LBM")
    ax.plot(r, u_ana, lw=1.8, ls="--", label="Analytical")
    ax.set_title(f"{title} | RMS = {rms(err):.3e}")
    ax.grid(True)
    ax.legend()

fig.supxlabel(r"$r$")
fig.supylabel(r"$u_\theta$")
plt.suptitle("Taylor–Couette: LBM vs Analytical")
plt.tight_layout(rect=[0, 0, 1, 0.95])

plt.savefig(out_dir / "u_theta_profiles_lbm_vs_analytical.png", dpi=300)
plt.close()

# ==================================================
# SANITY CHECK
# ==================================================
U_wall = Omega_i * Ri
print(f"Plots salvos em: {out_dir.resolve()}")
print(f"Expected wall velocity U_wall = {U_wall:.6g}")
print(f"LBM near inner wall (theta=0): u_0[0] = {u_0[0]:.6g}")
