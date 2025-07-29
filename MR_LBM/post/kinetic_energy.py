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

# Define simulation IDs and custom RGB colors
sim_ids = ["018", "019"]
labels = ["3rd","4th"]  # Consider using more descriptive labels
colors = [green, yellow, orange, red, purpple, blue]
linestyles = [(0, (1,1,1,2,6,2)), (0, (4,4)), (0, (1,2)), (0, (10,4)), (0, (6,2,1,2))]

# Create figure with better quality settings
plt.figure(figsize=(10, 6), dpi=150)

for sim_id, color, ls, label in zip(sim_ids, colors, linestyles, labels):
    # Construct file path for current simulation
    tke_path = f'LDC/{sim_id}/total_kinetic_energy.bin'
    
    try:
        # Read binary data
        with open(tke_path, 'rb') as f:
            data = np.frombuffer(f.read(), dtype=np.float32)
        
        # Extract time and kinetic energy
        t_star = data[::2]    # Even indices
        tke_sum = data[1::2]  # Odd indices
        
        t_star_clean = []
        tke_sum_clean = []
        
        maxValue = 0.0
        
        for t, tke in zip(t_star,tke_sum):
            if(t >= maxValue):
                maxValue = t
                t_star_clean.append(t)
                tke_sum_clean.append(tke)
                
        # Plot with RGB color
        plt.plot(t_star_clean, tke_sum_clean, 
                 color=color, 
                 linestyle=ls,
                 linewidth=2,
                 label=f'Order = {label}')
    
    except FileNotFoundError:
        print(f"Warning: File not found - {tke_path}")
        

# Configure plot with professional styling
plt.xlabel('Normalized Time (t*)', fontproperties=font_prop, fontweight='bold')
plt.ylabel('Total Kinetic Energy', fontproperties=font_prop, fontweight='bold')
plt.grid(alpha=0.2, linestyle='--')

# Create legend with custom font
legend = plt.legend(framealpha=0.9)
for text in legend.get_texts():
    text.set_fontproperties(font_prop)

# Add minor ticks and improve layout
plt.minorticks_on()
plt.tight_layout()

# Save as high-quality PDF (vector format)
plt.savefig('kinetic_energy_plot_high_order.pdf', format='pdf', bbox_inches='tight')

plt.show()