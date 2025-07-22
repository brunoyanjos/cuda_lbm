import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Create a font properties object for Times New Roman, falling back to DejaVu Serif
try:
    # First try to use Times New Roman
    font_prop = fm.FontProperties(family='Times', size=16)
    # Test if the font is available
    if not any(f.name == 'Times' for f in fm.fontManager.ttflist):
        raise ValueError("Times not found")
except:
    # Fallback to DejaVu Serif if Times is not available
    font_prop = fm.FontProperties(family='DejaVu Serif', size=16)

# Set global font properties
plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = font_prop.get_name()
plt.rcParams['mathtext.fontset'] = 'stix'  # For math text

# Caminho do arquivo
file_path = "LDC/012/velocity_probes.bin"

# Define o tipo de dado: use float32 se seu dfloat for float, float64 se for double
dtype = np.float32  # ou np.float32, se for o caso

# Parâmetros do arquivo
num_probes = 9
record_size = 1 + num_probes  # 1 para t_star + 9 sondas

# Leitura dos dados binários
data = np.fromfile(file_path, dtype=dtype)

# Reorganiza os dados por linha
data = data.reshape((-1, record_size))

# Separa tempo e valores das sondas
t_star = data[:, 0]
probes = data[:, 1:]

# Configurar figura com 2 subplots
plt.figure(figsize=(6, 6), dpi=300)

for i in range(num_probes):
    plt.plot(t_star, probes[:, i] / 0.0256, label=f"Sonda {i+1}")
    
plt.xlabel("t*")
plt.ylabel("$u_x$/$U_{lid}$")

# plt.legend()
plt.grid(True)
plt.tight_layout()

# Save as high-quality PDF (vector format)
plt.savefig('probes_100000.pdf', format='pdf', bbox_inches='tight')

plt.show()
    

