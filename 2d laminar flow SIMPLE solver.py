import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

# ==============================================================================
# 1. INPUTS: geometry, fluid properties, simulation parameters - edit to your liking or needs
# ==============================================================================

H = 0.1    #Height of pipe [m]

inlet_velocity = 0.01    # [m/s]

#material information 

rho = 980  # Density [kg/m^3]
mu_dynamic = 0.001   # Dynamic viscosity [kg/(m*s)]

#right now water^

# ==============================================================================
#Iteration parameters - increase for better precision
max_iter = 1000
#precision of convergence
tolerance = 1e-6
# ==============================================================================

# ==============================================================================
# 2. GEOMETRY CALCULATION: length of pipe for the velocity profile to fully develop and discretization
# ==============================================================================
# Reynolds number
Re = (rho * inlet_velocity * H) / mu_dynamic

#laminar flow check
if Re >= 2300:
    exit(f"ERROR: Reynolds number is {Re:.2f} Flow is not laminar! Please reduce inlet velocity or increase viscosity. Reynolds number must be < 2300 for this simulation to be valid.")
else:   print(f"Reynolds number is {Re:.2f} - Flow is laminar, proceeding with simulation.")

# Hydrodynamic Entrance Length (L_e) for a 2D channel - lengthe needed for the velocity profile to fully develop
L_e = 0.05 * Re * H

#The total pipe length to be L_e plus a 20% 
L = L_e * 1.2 

# For creeping flow - condition of minimal length 
if L < 5 * H:
    L = 5 * H

print(f"--- PHYSICS SETUP ---")
print(f"Height of pipe: {H:.2f}")
print(f"Entrance Length required for the velocity profile to develop: {L_e:.2f} m")
print(f"Setting Pipe Length (L) to: {L:.2f} m\n")

# Discretization resolution - grid
nx = 100   
ny = 40    

#cell size
dx = L / (nx - 1)
dy = H / (ny - 1)

#mesh
x = np.linspace(0, L, nx)
y = np.linspace(0, H, ny)
X, Y = np.meshgrid(x, y)


# Initial guess of velocity for all nodes same as inlet velocity 
u = np.ones((ny, nx)) * inlet_velocity
# Force the walls back to 0 m/s (No-slip condition)
u[0, :] = 0.0
u[-1, :] = 0.0
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))

#relaxation coefficients for pressure and velocity used in SIMPLE algorithm to "fight divergence"
if Re < 700 or mu_dynamic < 0.001:
    alpha_u = 0.3
    alpha_v = 0.4
    alpha_p = 0.1
else:
    alpha_u = 0.4
    alpha_v = 0.5
    alpha_p = 0.3   

# ==============================================================================
# 3. SOLVER FUNCTIONS
# ==============================================================================
# momentum x 
def solve_momentum_x(u, v, p, rho, mu, dx, dy, alpha_u, U_inlet):
    #setting up Ax = b, x = new velocity for each cell calculated based on pressure pressure field from previous step
    ny, nx = u.shape
    A = lil_matrix((nx*ny, nx*ny))
    b = np.zeros(nx*ny)

    #picking the cells by j and i through the grid
    for j in range(ny):
        for i in range(nx):
            row = i + j * nx
            
            if j == 0 or j == ny - 1: # Walls
                A[row, row] = 1.0
                b[row] = 0.0
            elif i == 0:              # Inlet
                A[row, row] = 1.0
                b[row] = U_inlet
            elif i == nx - 1:         # Outlet
                A[row, row] = 1.0
                A[row, row - 1] = -1.0
                b[row] = 0.0
            else:                     # Interior cells

                #convective mass fluxes for cells around centre cell
                Fe = rho * 0.5 * (u[j, i+1] + u[j, i]) * dy
                Fw = rho * 0.5 * (u[j, i-1] + u[j, i]) * dy
                Fn = rho * 0.5 * (v[j+1, i] + v[j, i]) * dx
                Fs = rho * 0.5 * (v[j-1, i] + v[j, i]) * dx

                #diffusive coeficients for cells around centre cell
                De = mu * dy / dx
                Dw = mu * dy / dx
                Dn = mu * dx / dy
                Ds = mu * dx / dy

                # sum of convective and diffusive terms
                ae = De + max(0, -Fe)
                aw = Dw + max(0, Fw)
                an = Dn + max(0, -Fn)
                as_ = Ds + max(0, Fs)
                
                #sum of cells around the centre cell - momentum 
                ap = ae + aw + an + as_ 

                #brakes on sharp gradients             
                ap_relaxed = ap / alpha_u 

                A[row, row] = ap_relaxed
                A[row, row + 1] = -ae
                A[row, row - 1] = -aw
                A[row, row + nx] = -an
                A[row, row - nx] = -as_

                #current pressure geradient in x dir
                pressure_grad = (p[j, i-1] - p[j, i+1]) * 0.5 * dy 
                #relaxed pressure changing term
                relaxation_source = (1 - alpha_u) * ap_relaxed * u[j, i]

                b[row] = pressure_grad + relaxation_source

    #solve system of matrices
    u_star_flat = spsolve(A.tocsr(), b)
    #new iterated velocity
    return u_star_flat.reshape((ny, nx)), A.diagonal() 

def solve_momentum_y(u, v, p, rho, mu, dx, dy, alpha_v):
    ny, nx = v.shape
    A = lil_matrix((nx*ny, nx*ny))
    b = np.zeros(nx*ny)
    for j in range(ny):
        for i in range(nx):
            row = i + j * nx
            if j == 0 or j == ny - 1:
                A[row, row] = 1.0
                b[row] = 0.0
            elif i == 0:
                A[row, row] = 1.0
                b[row] = 0.0
            elif i == nx - 1:
                A[row, row] = 1.0
                A[row, row - 1] = -1.0
                b[row] = 0.0
            else:
                Fe = rho * 0.5 * (u[j, i+1] + u[j, i]) * dy
                Fw = rho * 0.5 * (u[j, i-1] + u[j, i]) * dy
                Fn = rho * 0.5 * (v[j+1, i] + v[j, i]) * dx
                Fs = rho * 0.5 * (v[j-1, i] + v[j, i]) * dx
                De = mu * dy / dx; Dw = mu * dy / dx
                Dn = mu * dx / dy; Ds = mu * dx / dy
                ae = De + max(0, -Fe)
                aw = Dw + max(0, Fw)
                an = Dn + max(0, -Fn)
                as_ = Ds + max(0, Fs)
                ap = ae + aw + an + as_ 
                ap_relaxed = ap / alpha_v
                A[row, row] = ap_relaxed
                A[row, row + 1] = -ae
                A[row, row - 1] = -aw
                A[row, row + nx] = -an
                A[row, row - nx] = -as_
                pressure_grad = (p[j-1, i] - p[j+1, i]) * 0.5 * dx 
                relaxation_source = (1 - alpha_v) * ap_relaxed * v[j, i]
                b[row] = pressure_grad + relaxation_source
    v_star_flat = spsolve(A.tocsr(), b)
    return v_star_flat.reshape((ny, nx)), A.diagonal()

#balancing the forces acting on the fluid with mass flux across boundaries to get more accurate pressure field - fulfilling the continuity equation
def solve_pressure_correction(u_star, v_star, rho, dx, dy, ap_u, ap_v):
    #setting up solution
    ny, nx = u_star.shape
    A = lil_matrix((nx*ny, nx*ny))
    b = np.zeros(nx*ny)

    # moment of inertia in x and y direction for each cell from momentum equations
    ap_u_2d = ap_u.reshape((ny, nx))
    ap_v_2d = ap_v.reshape((ny, nx))

    #sorting through the cells by j and i
    for j in range(ny):
        for i in range(nx):
            row = i + j * nx
            
            if i == nx - 1: #inlet
                A[row, row] = 1.0
                b[row] = 0.0
            elif i == 0 or j == 0 or j == ny - 1: #walls, outlet    
                A[row, row] = 1.0
                if i == 0: A[row, row + 1] = -1.0
                elif j == 0: A[row, row + nx] = -1.0
                elif j == ny - 1: A[row, row - nx] = -1.0
                b[row] = 0.0
            else:
                #presssure coefficients based on momentum equations for each cell
                de = dy * dy / ap_u_2d[j, i+1] if i+1 < nx else 0
                dw = dy * dy / ap_u_2d[j, i-1] if i-1 >= 0 else 0
                dn = dx * dx / ap_v_2d[j+1, i] if j+1 < ny else 0
                ds = dx * dx / ap_v_2d[j-1, i] if j-1 >= 0 else 0
                
                #moment of inertia from momentum equations for each cell
                ae = rho * de; aw = rho * dw
                an = rho * dn; as_ = rho * ds
                ap = ae + aw + an + as_

                A[row, row] = ap
                A[row, row + 1] = -ae
                A[row, row - 1] = -aw
                A[row, row + nx] = -an
                A[row, row - nx] = -as_

                #mass flux in each cell 
                mass_flux_out = rho * (u_star[j, i+1] - u_star[j, i-1]) * 0.5 * dy + \
                                rho * (v_star[j+1, i] - v_star[j-1, i]) * 0.5 * dx
                b[row] = -mass_flux_out

    p_corr_flat = spsolve(A.tocsr(), b)
    return p_corr_flat.reshape((ny, nx))

#correctiong the velocity field based on the new pressure
def correct_u(u_star, p_corr, dx, dy, ap_u, alpha_p):
    #setup
    ny, nx = u_star.shape
    u_new = np.copy(u_star)
    ap_u_2d = ap_u.reshape((ny, nx))
    for j in range(1, ny-1):
        for i in range(1, nx-1):
            #pressure correction gradient in x direction
            d_p_dx = (p_corr[j, i-1] - p_corr[j, i+1]) * 0.5 / dx
            #inertia term 
            d_coef = dy * dx / ap_u_2d[j, i]
            #new velocity
            u_new[j, i] = u_star[j, i] + alpha_p * d_coef * d_p_dx
    return u_new

def correct_v(v_star, p_corr, dx, dy, ap_v, alpha_p):
    ny, nx = v_star.shape
    v_new = np.copy(v_star)
    ap_v_2d = ap_v.reshape((ny, nx))
    for j in range(1, ny-1):
        for i in range(1, nx-1):
            d_p_dy = (p_corr[j-1, i] - p_corr[j+1, i]) * 0.5 / dy
            d_coef = dx * dy / ap_v_2d[j, i]
            v_new[j, i] = v_star[j, i] + alpha_p * d_coef * d_p_dy
    return v_new

# ==============================================================================
# 4. MAIN EXECUTION LOOP
# ==============================================================================
print("Starting SIMPLE loop...")

#values from previous iteration
u_old = np.copy(u)
v_old = np.copy(v)
p_old = np.copy(p)

#plot 
plt.ion()
fig = plt.figure()
history_u, history_v, history_p = [], [], []

for it in range(max_iter):
    #momentum
    u_star, ap_u = solve_momentum_x(u_old, v_old, p_old, rho, mu_dynamic, dx, dy, alpha_u, inlet_velocity)
    v_star, ap_v = solve_momentum_y(u_old, v_old, p_old, rho, mu_dynamic, dx, dy, alpha_u)
    #pressure corrrection
    p_corr = solve_pressure_correction(u_star, v_star, rho, dx, dy, ap_u, ap_v)

    #old perssure and pressure correction
    p_new = p_old + alpha_p * p_corr
    #velocity correction
    u_new = correct_u(u_star, p_corr, dx, dy, ap_u, alpha_p)
    v_new = correct_v(v_star, p_corr, dx, dy, ap_v, alpha_p)

    #residuals 
    error_u = np.linalg.norm(u_new - u_old)
    error_v = np.linalg.norm(v_new - v_old)
    error_p = np.linalg.norm(p_new - p_old)

    history_u.append(error_u)
    history_v.append(error_v)
    history_p.append(error_p)

    #max error for corvengence criteria 
    max_error = max(error_u, error_v, error_p)   

    #every 10 iter 
    if it % 10 == 0:
        print(f"Iter {it:4d} | Res_u: {error_u:.2e} | Res_v: {error_v:.2e} | Res_p: {error_p:.2e}")
        plt.clf()
        plt.plot(history_u, label='Velocity U')
        plt.plot(history_v, label='Velocity V')
        plt.plot(history_p, label='Pressure')
        
        plt.yscale('log')
        plt.xlabel('Iteration')
        plt.ylabel('Residual')
        plt.title(f'Residual Monitor (Iteration {it})')
        plt.legend()
        plt.grid(True, which="both", linestyle="--", alpha=0.5)
        
        plt.pause(0.01)

    if max_error < tolerance:
        print(f"--> CONVERGED successfully at iteration {it}!")
        break
        
    u_old = np.copy(u_new)
    v_old = np.copy(v_new)
    p_old = np.copy(p_new)
    
plt.ioff()
plt.close('all')
# ==============================================================================
# 5. VISUALIZATION
# ==============================================================================

if L_e < L:
    
    #finding x near L_e     
    develop_idx = np.searchsorted(x, L_e)
    
    # extracting u values at x
    u_developed = u_new[:, develop_idx]
    
    print(f"Flow fully developed at node {develop_idx} (X = {x[develop_idx]:.2f} m)")

print("Plotting results...")

#velocity vector magnitude 
velocity_magnitude = np.sqrt(u_new**2 + v_new**2)

# 3 subplots 
fig2, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9))

# Plot 1: Velocity Field --------
cp1 = ax1.contourf(X, Y, velocity_magnitude, levels=20, cmap='viridis', alpha=0.9)
fig2.colorbar(cp1, ax=ax1, label='Velocity Magnitude (m/s)')

# number of printed profiles
num_stations = 5

# Endpoint=False so it stops before hitting the outlet boundary
station_indices = np.linspace(0, nx - 1, num_stations, endpoint=False, dtype=int)

# Scale of profile
profile_scale = (L / num_stations) * 0.6 / inlet_velocity 

# Draw a velocity profile at each station
for i in station_indices:
    x_station = x[i]
    u_profile = u_new[:, i]
    
    #curve
    x_curve = x_station + (u_profile * profile_scale)
    
    # baseline of profile (where velocity = 0)
    ax1.plot([x_station, x_station], [0, H], color='black', linestyle='--', linewidth=1, alpha=0.6)
    
    # solid velocity profile curve
    ax1.plot(x_curve, y, color='black', linewidth=2)
    
    # area shading
    ax1.fill_betweenx(y, x_station, x_curve, color='black', alpha=0.3)

ax1.set_title("Velocity Field with developing velocity profile")
ax1.set_xlabel("x [m]")
ax1.set_ylabel("y [m]")
ax1.set_xlim(0, L)
ax1.set_ylim(0, H)

# 2. Plot - pressure field ------
cp2 = ax2.contourf(X, Y, p_new, levels=50, cmap='coolwarm')
fig2.colorbar(cp2, ax=ax2, label='Pressure (Pa relative to outlet)')
ax2.set_title("Pressure Field")
ax2.set_xlabel("x [m]")
ax2.set_ylabel("y [m]")

# Plot 3: velocity profile at outlet compared to analytical solution for fully developed flow --------
u_outlet = u_new[:, -1]

#analytical solution for 3rd plot
u_max_analytical = 1.5 * inlet_velocity
u_analytical = u_max_analytical * (1 - ((y - (H/2)) / (H/2))**2)

ax3.plot(u_outlet, y, 'b-', linewidth=3, label='Numerical Solver Result')
ax3.plot(u_analytical, y, 'r--', linewidth=2, label='Analytical Exact Solution')
ax3.set_title("Fully Developed Velocity Profile")
ax3.set_xlabel("u(y) [m/s]")
ax3.set_ylabel("y [m]")
ax3.legend()
ax3.grid(True)

plt.tight_layout()
plt.show()