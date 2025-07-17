import numpy as np
import matplotlib.pyplot as plt

# Carregar dados do arquivo .dat
data = np.loadtxt('GRID/001/macr_0.dat', skiprows=1)  # Substitua pelo nome do seu arquivo
x = data[:, 0]
y = data[:, 1]
rho = data[:, 2]
ux = data[:, 3]
uy = data[:, 4]

# Identificar valores únicos para reconstruir a malha
unique_x = np.unique(x)
unique_y = np.unique(y)
nx = len(unique_x)
ny = len(unique_y)

# Reformatar dados em grades 2D
X = x.reshape((ny, nx))
Y = y.reshape((ny, nx))
Rho = rho.reshape((ny, nx))
Ux = ux.reshape((ny, nx))
Uy = uy.reshape((ny, nx))

# Calcular magnitude da velocidade
velocity_magnitude = np.sqrt(Ux**2 + Uy**2)

# Criar figuras
plt.figure(figsize=(15, 10))

# Plot 1: Densidade (rho)
plt.subplot(2, 2, 1)
plt.pcolormesh(X, Y, Rho, shading='auto', cmap='viridis')
plt.colorbar(label='Densidade (rho)')
plt.title('Densidade')
plt.xlabel('x')
plt.ylabel('y')

# Plot 2: Magnitude da Velocidade
plt.subplot(2, 2, 2)
plt.pcolormesh(X, Y, velocity_magnitude, shading='auto', cmap='jet')
plt.colorbar(label='Magnitude da Velocidade')
plt.title('Magnitude da Velocidade')
plt.xlabel('x')
plt.ylabel('y')

# Plot 3: Campo Vetorial (reduzido para clareza)
plt.subplot(2, 2, 3)
skip = 2  # Reduzir número de setas
plt.quiver(
    X[::skip, ::skip], 
    Y[::skip, ::skip], 
    Ux[::skip, ::skip], 
    Uy[::skip, ::skip],
    scale=0.3,  # Ajuste conforme necessário
    color='white'
)
plt.title('Campo Vetorial de Velocidade')
plt.xlabel('x')
plt.ylabel('y')

# Plot 4: Linhas de Corrente
plt.subplot(2, 2, 4)
plt.streamplot(
    X, Y, Ux, Uy, 
    density=2.0, 
    color='black', 
    linewidth=1
)
plt.title('Linhas de Corrente')
plt.xlabel('x')
plt.ylabel('y')

plt.tight_layout()
plt.savefig('resultados_cavidade.png', dpi=300)
plt.show()