import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --------------------------------------------------
# Configurações
# --------------------------------------------------
path = Path("ANNUL/001/kinetic_energy.bin")
dtype = np.float32  # ajuste se for float32
dt = 1.0  # passo de tempo (LBM step = 1 normalmente)

# --------------------------------------------------
# Leitura do arquivo binário
# --------------------------------------------------
ke = np.fromfile(path, dtype=dtype)

time = np.arange(len(ke)) * dt

# --------------------------------------------------
# Plot
# --------------------------------------------------
plt.figure(figsize=(8, 4))
plt.plot(time, ke, lw=1.5)
plt.xlabel("Time [LBM steps]")
plt.ylabel("Kinetic Energy")
plt.title("Kinetic Energy vs Time")
plt.grid(True)
plt.tight_layout()
plt.show()
