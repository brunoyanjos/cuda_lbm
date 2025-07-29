import pandas as pd
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

# Define custom colors
green = '#61BB46'
yellow = '#FDB827'
orange = '#F5821F'
red = '#E03A3E'
purpple = '#963D97'
blue = '#009DDC'

sim_ids = ['018', '019']
labels = ['3rd', '4th']  # Consider using more descriptive labels
linestyles = [(0, (4,4)), 'solid']

# Carregar dados de referência
benchmark_data_x = pd.read_csv('MR_LBM/post/benchmark/ghia_ux_dy.csv', header=None)
benchmark_data_y = pd.read_csv('MR_LBM/post/benchmark/ghia_uy_dx.csv', header=None)

# Extrair dados para perfis verticais (u_x vs y)
y_d = benchmark_data_x.iloc[1:, 0].astype(float).values
u_x_bench = benchmark_data_x.iloc[1:, 1:].astype(float).values.T
re_numbers = [f"Re= {int(value)}" for value in benchmark_data_x.iloc[0, 1:].values]

# Extrair dados para perfis horizontais (u_y vs x)
x_d = benchmark_data_y.iloc[1:, 0].astype(float).values
u_y_bench = benchmark_data_y.iloc[1:, 1:].astype(float).values.T

benchmark_index = 6  # Re=10000

# Configurar figura com 2 subplots
plt.figure(figsize=(6, 6), dpi=300)

for sim_id, ls, label in zip(sim_ids, linestyles, labels):
    # Ler parâmetros da simulação
    info_path = f"LDC/{sim_id}/info.txt"
    try:
        with open(info_path, 'r') as f:
            params = {line.split(':')[0].strip(): line.split(':')[1].strip()
                      for line in f if ':' in line}
        nx = int(params.get('NX', 100))
        ny = int(params.get('NY', 100))
        umax = float(params.get('Umax', 1.0))
    except Exception as e:
        print(f"Erro em {sim_id}: {e}")
        continue
    
    # Construct file path for current simulation
    file_x_path = f'LDC/{sim_id}/velocity_avg_x.bin'
    file_y_path = f'LDC/{sim_id}/velocity_avg_y.bin'
    
    try:
        ux = np.fromfile(file_x_path, dtype=np.float32)
        uy = np.fromfile(file_y_path, dtype=np.float32)
            
        pos = np.linspace(0, 1, nx)
                            
        # # Plot with RGB color
        plt.plot(ux  / (2 * umax), pos - 1/2, 
                 color=blue, 
                 linestyle=ls,
                 linewidth=2,
                )
        
        plt.plot(pos - 1/2, uy / (2 * umax),
                 color=red, 
                 linestyle=ls,
                 linewidth=2,
                )
        
        dummy_handle, = plt.plot([], [], linestyle= ls, color='black', label=f"Order: {label}")
    
    except FileNotFoundError:
        print(f"Warning: File not found - {sim_id}") 

#Add benchmark plot and its handle to the lists
benchmark_handle_ux, = plt.plot(u_x_bench[benchmark_index] / 2, y_d - 0.5, 'x', color='black')

# === PERFIL HORIZONTAL (u_y vs x) ===
# No need to add benchmark_handle_uy as it's the same style as benchmark_handle_ux
plt.plot(x_d - 0.5, u_y_bench[benchmark_index] / 2, 'x', color='black')

dummy_handle, = plt.plot([], [], 'x', linestyle= '',  color='black', label=f"Benchmark")

#Configure plot with professional styling
plt.grid(alpha=0.2, linestyle='--')

# Create legend with custom font
legend = plt.legend(framealpha=0.9)
for text in legend.get_texts():
    text.set_fontproperties(font_prop)

# Add minor ticks and improve layout
plt.minorticks_on()
plt.tight_layout()

# Save as high-quality PDF (vector format)
plt.savefig('average_velocity_high_order.pdf', format='pdf', bbox_inches='tight')

plt.show()


