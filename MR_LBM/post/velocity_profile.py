import numpy as np
import matplotlib.pyplot as plt
import os

# =============================
# CONFIGURAÇÃO DO CAMINHO
# =============================

PATH = "ANNUL/001/"  # << ALTERE PARA O SEU DIRETÓRIO >>
UX_FILE = os.path.join(PATH, "ux_dy_average.bin")

DTYPE = np.float32  # mude para np.float64 se dfloat = double
NUM_POINTS = 9  # 3 × 3 probes

# =============================
# FUNÇÃO PARA LER O ARQUIVO
# =============================


def load_probe_series(filepath, dtype):
    """Carrega arquivo binário escrito como blocos de 9 floats."""
    data = np.fromfile(filepath, dtype=dtype)

    return data


print("Carregando arquivos...")

ux = load_probe_series(UX_FILE, DTYPE)

y_values = np.arange(0, 512, 1)

plt.plot(ux, y_values)

plt.show()
