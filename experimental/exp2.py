import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os

# ==========================================
# 1. ENGINEERING CONSTANTS
# ==========================================
OUTPUT_DIR = "sim_results_shapes"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Geometry
LENGTH = 20.0           # 20m Length to allow full flow recovery
R_PIPE_BASE = 0.5       
G = 9.81
PIPE_SLOPE = np.radians(30) # 30 deg slope

# Material (Savage-Hutter)
RHO_IN = 833.0          
PHI = np.radians(35)    
DELTA = np.radians(25)  
EPSILON = 1.0           

V_IN = 10.0              # Inlet Velocity

# ==========================================
# 2. PHYSICS ENGINE (Savage-Hutter Energy)
# ==========================================
def solve_scenario(angle_deg, lag_m, shape_type):
    # --- A. Geometry Setup ---
    SPLIT_START = 2.0
    SPLIT_LEN = 1.0
    split_max_r = np.tan(np.radians(angle_deg)) * (SPLIT_LEN / 2.0)
    
    # CONSTRAINT: Splitter cannot block > 85% of pipe
    if split_max_r >= (0.85 * R_PIPE_BASE):
        return 999.0, None, None, None 

    # Constriction Setup
    CONST_START = (SPLIT_START + SPLIT_LEN) + lag_m
    CONST_LEN = 1.0
    CONST_SQUEEZE = 0.15

    def get_geo(x):
        # 1. Splitter (Inner)
        r_in = 0.0
        if SPLIT_START <= x <= (SPLIT_START + SPLIT_LEN):
            lx = (x - SPLIT_START) / SPLIT_LEN
            if lx < 0.5: r_in = split_max_r * (lx * 2)
            else:        r_in = split_max_r * (2 * (1 - lx))
            
        # 2. Outer Wall (The Shape Variable)
        r_out = R_PIPE_BASE
        
        if CONST_START <= x <= (CONST_START + CONST_LEN):
            lx = (x - CONST_START) / CONST_LEN
            squeeze = 0.0
            
            # --- SHAPE LOGIC ---
            if shape_type == "Sine":
                squeeze = CONST_SQUEEZE * np.sin(lx * np.pi)
                
            elif shape_type == "Plateau":
                # Ramp 20%, Hold 60%, Release 20%
                if lx < 0.2:   squeeze = CONST_SQUEEZE * (lx / 0.2)
                elif lx > 0.8: squeeze = CONST_SQUEEZE * ((1.0 - lx) / 0.2)
                else:          squeeze = CONST_SQUEEZE
                
            elif shape_type == "SharkFin":
                # Fast Ramp Down (20%), Slow Release (80%)
                if lx < 0.2: squeeze = CONST_SQUEEZE * (lx / 0.2)
                else:        squeeze = CONST_SQUEEZE * ((1.0 - lx) / 0.8)

            elif shape_type == "Gaussian":
                # Bell Curve
                base = np.exp(-((lx - 0.5)**2) / (2 * 0.15**2))
                off = np.exp(-((0.0 - 0.5)**2) / (2 * 0.15**2))
                squeeze = CONST_SQUEEZE * (base - off) / (1 - off)
                
            r_out = R_PIPE_BASE - squeeze
            
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

    # --- C. Coefficients ---
    rad = np.sqrt(1 - (1 + np.tan(DELTA)**2) * np.cos(PHI)**2)
    k_base = 2 / np.cos(PHI)**2
    K_PASS = k_base * (1 + rad) - 1
    K_ACT = k_base * (1 - rad) - 1

    # --- D. Energy Integration ---
    u = np.zeros_like(x_pts)
    u[0] = V_IN
    
    for i in range(len(x_pts)-1):
        k_eff = K_PASS if dh_dx[i] > 0 else K_ACT
        
        # Accelerations (m/s^2)
        a_grav = G * np.sin(PIPE_SLOPE)
        a_fric = G * np.cos(PIPE_SLOPE) * np.tan(DELTA)
        a_press = EPSILON * k_eff * np.cos(PIPE_SLOPE) * dh_dx[i] * G
        
        a_net = a_grav - a_fric - a_press
        
        # Energy Step
        u_sq = u[i]**2 + 2 * a_net * dx
        
        if u_sq <= 0.1: 
            u[i+1] = 0.1
        else:
            u[i+1] = np.sqrt(u_sq)
            
    return u[-1], x_pts, u, get_geo

# ==========================================
# 3. VISUALIZATION
# ==========================================
def save_simulation_image(shape, angle, lag, v_out, x, u, geo):
    r_in = [geo(xi)[0] for xi in x]
    r_out = [geo(xi)[1] for xi in x]

    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Walls
    ax.plot(x, r_out, 'k-', lw=3, label="Outer Wall")
    ax.plot(x, [-r for r in r_out], 'k-', lw=3)
    ax.fill_between(x, r_in, [-r for r in r_in], color='silver', label='Splitter')
    
    # Streamlines
    n_lines = 30
    for k in range(n_lines + 1):
        pct = k/n_lines
        r_line = np.sqrt(pct * (np.array(r_out)**2 - np.array(r_in)**2) + np.array(r_in)**2)
        
        points = np.array([x, r_line]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        # Fixed Norm 0-15 m/s
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
    
    ax.set_title(f"Shape: {shape} | Angle: {angle:.1f}° | Lag: {lag:.1f}m | Exit V: {v_out:.2f} m/s")
    ax.set_ylim(-0.8, 0.8)
    ax.set_xlim(0, LENGTH)
    ax.set_xlabel("Axial Distance (m)")
    ax.grid(True, alpha=0.3)
    
    filename = f"{OUTPUT_DIR}/{shape}_A{angle:.1f}_L{lag:.1f}.png"
    plt.savefig(filename)
    plt.close(fig)
    return filename

# ==========================================
# 4. MAIN OPTIMIZATION LOOP
# ==========================================
print(f"{'Shape':<10} | {'Angle':<8} | {'Lag':<8} | {'Exit Vel':<12} | {'Result'}")
print("-" * 60)

shapes = ["Sine", "Plateau", "SharkFin", "Gaussian"]
angles = [0,2.0,7.0,8.0,9.0,15.0,20.0,24.0,30.0,35.0,44.0]        # Test small to medium angles
lags = [0.2,0.3,2.0,3.0,4.0]          # Test mid to far range

global_best_v = 999.0
global_best_cfg = None

for shape in shapes:
    shape_best_v = 999.0
    
    for ang in angles:
        for lag in lags:
            v_out, x, u, geo = solve_scenario(ang, lag, shape)
            
            if v_out == 999.0: continue
            
            # Save Image
            save_simulation_image(shape, ang, lag, v_out, x, u, geo)
            
            print(f"{shape:<10} | {ang:<8.1f} | {lag:<8.1f} | {v_out:<12.2f} | Saved")
            
            # Track Bests (Must be > 0.5 to not be a clog)
            if v_out < global_best_v and v_out > 0.5:
                global_best_v = v_out
                global_best_cfg = (shape, ang, lag)
            
            if v_out < shape_best_v and v_out > 0.5:
                shape_best_v = v_out

    print(f"   >>> Best for {shape}: {shape_best_v:.2f} m/s")
    print("-" * 60)

print("OPTIMIZATION COMPLETE.")
if global_best_cfg:
    print(f"Global Best Design: {global_best_cfg[0]} Shape")
    print(f"Configuration: Angle {global_best_cfg[1]} deg, Lag {global_best_cfg[2]} m")
    print(f"Final Velocity: {global_best_v:.2f} m/s")
else:
    print("No valid designs found.")
