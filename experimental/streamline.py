import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os

# ==========================================
# 1. ENGINEERING CONSTANTS
# ==========================================
OUTPUT_DIR = "sim_results_all"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Geometry - 20.0m Length to ensure we see the full exit behavior
LENGTH = 50.0           
R_PIPE_BASE = 0.5       
G = 9.81
PIPE_SLOPE = np.radians(30) # 30 deg slope

# Material Properties (Savage-Hutter)
RHO_IN = 833.0          
PHI = np.radians(35)    
DELTA = np.radians(25)  
EPSILON = 1.0           

V_IN = 5.0              # Inlet Velocity

# ==========================================
# 2. PHYSICS ENGINE (Savage-Hutter Energy)
# ==========================================
def solve_scenario(angle_deg, lag_m):
    # --- A. Geometry Setup ---
    SPLIT_START = 2.0
    SPLIT_LEN = 10.0
    
    # Calculate Max Radius of Splitter based on Angle
    split_max_r = np.tan(np.radians(angle_deg)) * (SPLIT_LEN / 2.0)
    
    # CONSTRAINT: Splitter cannot block > 85% of pipe
    if split_max_r >= (0.85 * R_PIPE_BASE):
        return 999.0, None, None, None # INVALID

    # Constriction Setup
    CONST_START = (SPLIT_START + SPLIT_LEN) + lag_m
    CONST_LEN = 3.0
    CONST_SQUEEZE = 0.15 

    def get_geo(x):
        # 1. Splitter (Inner)
        r_in = 0.0
        if SPLIT_START <= x <= (SPLIT_START + SPLIT_LEN):
            lx = (x - SPLIT_START) / SPLIT_LEN
            if lx < 0.5: r_in = split_max_r * (lx * 2)
            else:        r_in = split_max_r * (2 * (1 - lx))
            
        # 2. Outer Wall
    
        r_out = R_PIPE_BASE
        if CONST_START <= x <= (CONST_START + CONST_LEN):
            lx = (x - CONST_START) / CONST_LEN
            if lx < 0.2:
                squeeze = CONST_SQUEEZE * (lx / 0.2)
            elif lx > 0.8:
                squeeze = CONST_SQUEEZE * ((1.0 - lx) / 0.2)
            else:
                squeeze = CONST_SQUEEZE
                r_out = R_PIPE_BASE - squeeze
        if CONST_START <= x <= (CONST_START + CONST_LEN):
            lx = (x - CONST_START) / CONST_LEN
            r_out = R_PIPE_BASE - (CONST_SQUEEZE * lx)
            
        return r_in, r_out

    # --- B. Virtual Depth Profile ---
    x_pts = np.linspace(0, LENGTH, 800)
    dx = x_pts[1] - x_pts[0]
    
    areas = []
    for x in x_pts:
        ri, ro = get_geo(x)
        areas.append(np.pi*(ro**2 - ri**2))
    areas = np.array(areas)
    
    h_virtual = areas[0] / areas 
    dh_dx = np.gradient(h_virtual, dx)

    # --- C. Savage-Hutter Coefficients ---
    rad = np.sqrt(1 - (1 + np.tan(DELTA)**2) * np.cos(PHI)**2)
    k_base = 2 / np.cos(PHI)**2
    K_PASS = k_base * (1 + rad) - 1
    K_ACT = k_base * (1 - rad) - 1

    # --- D. Integration ---
    u = np.zeros_like(x_pts)
    u[0] = V_IN
    
    for i in range(len(x_pts)-1):
        k_eff = K_PASS if dh_dx[i] > 0 else K_ACT
        
        a_grav = G * np.sin(PIPE_SLOPE)
        a_fric = G * np.cos(PIPE_SLOPE) * np.tan(DELTA)
        a_press = EPSILON * k_eff * np.cos(PIPE_SLOPE) * dh_dx[i] * G
        
        a_net = a_grav - a_fric - a_press
        
        u_sq = u[i]**2 + 2 * a_net * dx
        
        if u_sq <= 0.1: 
            u[i+1] = 0.1
        else:
            u[i+1] = np.sqrt(u_sq)
            
    return u[-1], x_pts, u, get_geo

# ==========================================
# 3. PLOTTING FUNCTION
# ==========================================
def save_simulation_image(angle, lag, v_out, x, u, geo):
    r_in = [geo(xi)[0] for xi in x]
    r_out = [geo(xi)[1] for xi in x]

    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Walls
    ax.plot(x, r_out, 'k-', lw=3, label="Outer Wall")
    ax.plot(x, [-r for r in r_out], 'k-', lw=3)
    ax.fill_between(x, r_in, [-r for r in r_in], color='silver', label='Splitter')
    
    # Streamlines
    n_lines = 10
    for k in range(n_lines + 1):
        pct = k/n_lines
        r_line = np.sqrt(pct * (np.array(r_out)**2 - np.array(r_in)**2) + np.array(r_in)**2)
        
        points = np.array([x, r_line]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        # Color Map (Fixed scale 0 to 15 m/s)
        norm = plt.Normalize(0, 15)
        lc = LineCollection(segments, cmap='jet_r', norm=norm)
        lc.set_array(u)
        lc.set_linewidth(1.5)
        ax.add_collection(lc)
        
        if k > 0:
            points_b = np.array([x, -r_line]).T.reshape(-1, 1, 2)
            segments_b = np.concatenate([points_b[:-1], points_b[1:]], axis=1)
            lc_b = LineCollection(segments_b, cmap='jet_r', norm=norm)
            lc_b.set_array(u)
            lc_b.set_linewidth(1.5)
            ax.add_collection(lc_b)

    cbar = plt.colorbar(lc, ax=ax)
    cbar.set_label('Velocity (m/s)', rotation=270, labelpad=15)
    
    ax.set_title(f"Config: Angle {angle:.1f}° | Lag {lag:.1f}m | Exit V: {v_out:.2f} m/s")
    ax.set_ylim(-0.8, 0.8)
    ax.set_xlim(0, LENGTH)
    ax.set_xlabel("Axial Distance (m)")
    ax.grid(True, alpha=0.3)
    
    # Unique Filename for EVERY iteration
    filename = f"{OUTPUT_DIR}/Design_A{angle:.1f}_L{lag:.1f}.png"
    plt.savefig(filename)
    plt.close(fig) # Close memory to prevent leaks
    return filename

# ==========================================
# 4. MAIN OPTIMIZATION LOOP
# ==========================================
print(f"{'Angle':<10} | {'Lag':<10} | {'Exit Vel':<15} | {'Result'}")
print("-" * 55)

angles = np.linspace(1, 6, 5)    # 5.0, 7.5, 10.0
lags = np.linspace(5.0, 9.0, 5)  # -1.0 to 5.0

best_v = 999.0
best_cfg = None

for ang in angles:
    for lag in lags:
        # Run Physics
        v_out, x, u, geo = solve_scenario(ang, lag)
        
        # Skip Invalid
        if v_out == 999.0:
            continue
            
        # *** SAVE IMAGE HERE (Inside Loop) ***
        save_simulation_image(ang, lag, v_out, x, u, geo)
        
        print(f"{ang:<10.1f} | {lag:<10.1f} | {v_out:<15.2f} | Image Saved")
        
        # Track Best
        if v_out < best_v and v_out > 0.5:
            best_v = v_out
            best_cfg = (ang, lag)

print("-" * 55)
if best_cfg:
    print(f"OPTIMIZATION COMPLETE.")
    print(f"Best Configuration: Angle {best_cfg[0]:.1f} deg, Lag {best_cfg[1]:.1f} m")
    print(f"Lowest Velocity: {best_v:.2f} m/s")
    print(f"All images saved in folder: {OUTPUT_DIR}/")
else:
    print("No valid design found.")
