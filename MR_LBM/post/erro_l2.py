import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from my_functions import tke_interpolation, find_idx
from scipy.stats import linregress

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

# Create figure with better quality settings
plt.figure(figsize=(6, 6), dpi=300)
    
base_sim = '005'

tke_path = f'LDC/{base_sim}/total_kinetic_energy.bin'
    
try:
    # Read binary data
    with open(tke_path, 'rb') as f:
        data = np.frombuffer(f.read(), dtype=np.float32)
        
    # Extract time and kinetic energy
    t_star_b = data[::2]   # Even indices
    tke_b = data[1::2]     # Odd indices
    
except FileNotFoundError:
    print(f"Warning: File not found - {tke_path}")
    
id_sims = ['001', '002', '003', '004']
grid_size = [128, 256, 512, 1024]

square_sum = sum(x**2 for x in tke_b)

l2_results = []

for id_sim in id_sims:
    tke_path = f'LDC/{id_sim}/total_kinetic_energy.bin'
    
    try:
        # Read binary data
        with open(tke_path, 'rb') as f:
            data = np.frombuffer(f.read(), dtype=np.float32)
        
        # Extract time and kinetic energy
        t_star = data[::2]   # Even indices
        tke = data[1::2]     # Odd indices
        
        l2 = 0.0
                  
        for t_ref, tke_ref in zip(t_star_b, tke_b):
            matches = find_idx(t_star, t_ref)
            
            if len(matches) == 1:
                idx = matches[0]
                
                l2 += (tke[idx] - tke_ref) ** 2
            elif len(matches) == 2:
                idx0 = matches[0]
                idx1 = matches[1]
                
                tke_int = tke_interpolation(t_star[idx0], tke[idx0], t_star[idx1], tke[idx1], t_ref)
                
                l2 += (tke_int - tke_ref) ** 2
        
        l2 = (l2 / square_sum) ** (0.5)
        
        l2_results.append(l2)
                          
    except FileNotFoundError:
        print(f"Warning: File not found - {tke_path}")
        
l2_results_clean = [float(x) for x in l2_results]

grid_log = np.log(grid_size)
l2_log = np.log(l2_results_clean)

slope_log, intercept_log, r_value, p_value, std_err = linregress(grid_log, l2_log)

b_estimado = slope_log
a_estimado = np.exp(intercept_log) # np.exp é a função exponencial (e^x)

y_regressao = a_estimado * (grid_size**b_estimado)

plt.loglog(grid_size, l2_results_clean, linestyle='', marker='o', color= 'black', label='Error')
plt.loglog(grid_size, y_regressao, linestyle='--', color= '#E03A3E', label='Regression')

# Configure plot with professional styling
plt.xlabel('Grid Size', fontproperties=font_prop, fontweight='bold')
plt.ylabel('$L_2$', fontproperties=font_prop, fontweight='bold')
plt.grid(True, which='both', alpha=0.2, linestyle='--')

plt.tick_params(axis='x', which='both', length=4.5, width=1.2, labelsize=12)
plt.tick_params(axis='y', which='both', length=4.5, width=1.2, labelsize=12)

# Create legend with custom font
legend = plt.legend(framealpha=0.9)
for text in legend.get_texts():
    text.set_fontproperties(font_prop)

# Add minor ticks and improve layout
plt.minorticks_on()
plt.tight_layout()

# Save as high-quality PDF (vector format)
plt.savefig('l2_error_Re_10000.pdf', format='pdf', bbox_inches='tight')
print(l2_results_clean)
print(b_estimado)

plt.show()      
