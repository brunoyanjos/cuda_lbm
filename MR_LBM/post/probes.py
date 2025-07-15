import numpy as np
import matplotlib.pyplot as plt

# Caminho do arquivo
file_path = "LDC/101/velocity_probes.bin"

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

# Plotagem das 9 sondas
plt.figure(figsize=(10, 6))
for i in range(num_probes):
    plt.plot(t_star, probes[:, i] / 0.0256, label=f"Sonda {i+1}")
    
plt.xlabel("t*")
plt.ylabel("Valor das sondas")
plt.title("Evolução temporal das sondas de velocidade")
# plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
    

