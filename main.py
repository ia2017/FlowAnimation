# --- Fluid simulation in a defined 2-D quadrilateral grid ----
# --- Imran bin Ahmad Azhar

from FluidDriver import Fluid

# --- Input values ----

Re=1000 # Reynolds number
dt=0.1 # Timestep
T=5.0 # Final Time
Lx=1  # Dimensions
Ly=1
Nx=50   # Number of grid points
Ny=50   # Number of grid points

dx=Lx/(Nx-1)
dy = Ly / (Ny - 1)
# Check if CFL condition is satisfied
try:
    dt_check=Re*dx*dy/4
except dt>=dt_check:
    print("CFL condition not met")


# --- Initiating Fluid class ----

sim=Fluid(Re,dt,T,Lx,Ly,Nx,Ny)
sim.SetGrid()
sim.BoundaryConditions()
sim.Integrate()
sim.PrintGrid()
#sim.Contour()
sim.Test()




