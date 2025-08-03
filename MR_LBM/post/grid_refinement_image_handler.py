import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ===== CONFIGURAÇÕES =====
run = '001'
max_tempo = 50          # Número máximo de frames (tempo)
salvar_gif = True
gif_nome = f"comparacao_{run}.gif"
pasta = os.path.join("GRID", run)

# ===== FUNÇÃO PARA LER ARQUIVOS .dat =====
def ler_dados(prefixo, max_tempo):
    dados = []
    tempos = []
    
    for i in range(max_tempo + 1):
        nome_arquivo = f"{prefixo}_{i}.dat"
        caminho = os.path.join(pasta, nome_arquivo)
        
        if not os.path.exists(caminho):
            break
        
        matriz = np.loadtxt(caminho)
        dados.append(matriz)
        tempos.append(i)
    
    return np.array(dados), tempos

# ===== LEITURA DE DADOS COARSE E FINE =====
dados_coarse, tempos = ler_dados("coarse_macr", max_tempo)
dados_fine, _        = ler_dados("fine_macr", max_tempo)

# Verificação
assert len(dados_coarse) == len(dados_fine), "Nº de frames diferentes entre coarse e fine!"

# ===== PLOTAGEM E ANIMAÇÃO =====
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

vmin = min(dados_coarse.min(), dados_fine.min())
vmax = max(dados_coarse.max(), dados_fine.max())

im1 = ax1.imshow(dados_coarse[0], cmap="seismic", origin="lower", vmin=vmin, vmax=vmax)
ax1.set_title("Coarse Grid")
plt.colorbar(im1, ax=ax1)

im2 = ax2.imshow(dados_fine[0], cmap="seismic", origin="lower", vmin=vmin, vmax=vmax)
ax2.set_title("Fine Grid")
plt.colorbar(im2, ax=ax2)

title = fig.suptitle(f"Tempo: {tempos[0]}")

def update(frame):
    im1.set_data(dados_coarse[frame])
    im2.set_data(dados_fine[frame])
    title.set_text(f"Tempo: {tempos[frame]}")
    return [im1, im2, title]

ani = animation.FuncAnimation(fig, update, frames=len(dados_coarse), blit=True)

plt.tight_layout()
plt.show()

# ===== SALVAR GIF OPCIONAL =====
if salvar_gif:
    ani.save(gif_nome, writer="pillow", fps=5)
    print(f"GIF salvo como: {gif_nome}")
