import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as fm

# Create a font properties object for Times New Roman, falling back to DejaVu Serif
try:
    # First try to use Times New Roman
    font_prop = fm.FontProperties(family="Times", size=16)
    # Test if the font is available
    if not any(f.name == "Times" for f in fm.fontManager.ttflist):
        raise ValueError("Times not found")
except:
    # Fallback to DejaVu Serif if Times is not available
    font_prop = fm.FontProperties(family="DejaVu Serif", size=16)

# Lista de IDs de simulação
sim_ids = ["001"]
labels = ["IRBC"]

# Define a color for the simulation lines
green = "#61BB46"
yellow = "#FDB827"
orange = "#F5821F"
red = "#E03A3E"
purpple = "#963D97"
blue = "#009DDC"

linestyles = [(0, (4, 4)), "solid", (0, (1, 2)), (0, (10, 4)), "solid"]
colors_x = [green, orange, red, purpple, blue]

# Carregar dados de referência
# benchmark_data_x = pd.read_csv('MR_LBM/post/benchmark/ghia_ux_dy.csv', header=None)
# benchmark_data_y = pd.read_csv('MR_LBM/post/benchmark/ghia_uy_dx.csv', header=None)

# Extrair dados para perfis verticais (u_x vs y)
# y_d = benchmark_data_x.iloc[1:, 0].astype(float).values
# u_x_bench = benchmark_data_x.iloc[1:, 1:].astype(float).values.T
# re_numbers = [f"Re= {int(value)}" for value in benchmark_data_x.iloc[0, 1:].values]

# Extrair dados para perfis horizontais (u_y vs x)
# x_d = benchmark_data_y.iloc[1:, 0].astype(float).values
# u_y_bench = benchmark_data_y.iloc[1:, 1:].astype(float).values.T

# Configurar figura com 2 subplots
plt.figure(figsize=(6, 6), dpi=150)

# Índice para número de Reynolds desejado
benchmark_index = 6  # Re=10000

# Lists to store handles and labels for the legend
legend_handles = []
legend_labels = []

# Processar cada simulação
for sim_id, line, color_x, label in zip(sim_ids, linestyles, colors_x, labels):
    # Ler parâmetros da simulação
    info_path = f"LDC/{sim_id}/info.txt"
    try:
        with open(info_path, "r") as f:
            params = {
                line.split(":")[0].strip(): line.split(":")[1].strip()
                for line in f
                if ":" in line
            }
        nx = int(params.get("NX", 100))
        ny = int(params.get("NY", 100))
        umax = float(params.get("Umax", 1.0))
    except Exception as e:
        print(f"Erro em {sim_id}: {e}")
        continue

    # Ler velocidades
    try:
        ux_sim = np.fromfile(f"LDC/{sim_id}/velocity_x.bin", dtype=np.float32)
        uy_sim = np.fromfile(f"LDC/{sim_id}/velocity_y.bin", dtype=np.float32)
    except Exception as e:
        print(f"Erro ao ler dados de {sim_id}: {e}")
        continue

    # Gerar coordenadas normalizadas
    y_sim = np.linspace(0, 1, ny)
    x_sim = np.linspace(0, 1, nx)

    # Adicionar aos gráficos
    # Plot u_x, but don't add label directly here
    plt.plot(ux_sim / (2 * umax), y_sim - 1 / 2, linestyle=line, color=blue)
    # Plot u_y, but don't add label directly here
    plt.plot(x_sim - 1 / 2, uy_sim / (2 * umax), linestyle=line, color=red)

    # Create a single dummy plot entry for the legend
    # This plot won't be visible, but its handle will have the correct linestyle and color
    (dummy_handle,) = plt.plot(
        [], [], linestyle=line, color="black", label=f"BC type: {label}"
    )
    legend_handles.append(dummy_handle)
    legend_labels.append(f"BC type: {label}")

# === PERFIL VERTICAL (u_x vs y) ===
# Add benchmark plot and its handle to the lists
# (benchmark_handle_ux,) = plt.plot(
#     u_x_bench[benchmark_index] / 2, y_d - 0.5, "x", color="black"
# )
(dummy_handle,) = plt.plot([], [], "x", linestyle="", color="black", label=f"Benchmark")

# === PERFIL HORIZONTAL (u_y vs x) ===
# No need to add benchmark_handle_uy as it's the same style as benchmark_handle_ux
# plt.plot(x_d - 0.5, u_y_bench[benchmark_index] / 2, "x", color="black")

# Configurar gráfico do perfil horizontal
plt.grid(True, which="both", alpha=0.2, linestyle="--")

plt.tick_params(axis="x", which="both", length=4.5, width=1.2, labelsize=12)
plt.tick_params(axis="y", which="both", length=4.5, width=1.2, labelsize=12)

# Use the collected handles and labels to create the legend
# Create legend with custom font
legend = plt.legend(framealpha=0.9)
for text in legend.get_texts():
    text.set_fontproperties(font_prop)

plt.tight_layout()
plt.savefig(
    "comparison_velocity_profiles_RBC.pdf", dpi=300, format="pdf", bbox_inches="tight"
)
plt.show()
