# 2D Laminar Pipe Flow CFD Solver (SIMPLE Algorithm)

This repository contains a python 2D Computational Fluid Dynamics (CFD) solver. It simulates the development of a laminar flow velocity profile inside a pipe using the **SIMPLE** (Semi-Implicit Method for Pressure Linked Equations) algorithm. 

## Functions

* **Customisable geometry, fluid medium and simulation precision:** In the first section of the code you can alter flow domain properties to your liking.
* **Laminar flow check:** Code checks to ensure the flow remains laminar (Reynolds number must be $< 2300$).
* **Automatic Length Setup:** Calculates the hydrodynamic entrance length needed for the velocity profile to fully develop based on the Reynolds number.
* **Residual Monitor:** A real-time plot tracks velocity and pressure errors to visually display the convergence of solution.
* **Visualization:** Generation of three subplots showing the velocity field with developing profiles, the pressure field, and a comparison of code numerical solution against the exact analytical solution.

---

## Prerequisites

To run this simulation, you will need Python installed along with the following libraries:

* `numpy` (Array operations and math)
* `matplotlib` (Visualizations and live plotting)
* `scipy` (Sparse matrices and linear algebra solvers)

---

## How to Use (Inputs)

The physics setup is conveniently located at the very top of the script in the **INPUTS** section. You can easily edit these to your liking or needs:

* **Geometry:** Set the height of the pipe (`H`). The length is handled automatically.
* **Fluid Properties:** Define the density (`rho`) and dynamic viscosity (`mu_dynamic`). Right now, the default material information is set to simulate water.
* **Boundary Conditions:** Set the initial `inlet_velocity`. 
* **Iteration Parameters:** You can increase `max_iter` for highly complex runs or adjust the `tolerance` ($1 \times 10^{-6}$ by default) for a stricter precision of convergence.

---

## How It Works (Under the Hood)

This script solves the 2D Navier-Stokes equations by discretizing the pipe into a grid ($100 \times 40$ cells by default) and iteratively solving the fluid flow matrices. Here is the step-by-step logic:

### 1. Geometry Calculation and Discretization
Before the main loop starts, the script calculates the Reynolds number.

$$Re = \frac{\rho \cdot U_{inlet} \cdot H}{\mu}$$

If the flow is safely laminar, it calculates the Hydrodynamic Entrance Length ($L_e$). This is the exact length needed for the velocity profile to fully develop. To be safe, the total pipe length is set to $L_e$ plus a 20% margin. The domain is then split into a mesh grid of cells. 

### 2. The SIMPLE Solver Loop
The algorithm relies on three main solver functions that repeat until the residual error hits the specified tolerance:

* **Solving Momentum ($x$ and $y$):** For each cell, we set up $Ax = b$ matrices based on the convective and diffusive fluxes from surrounding cells. This gives us a "star" (guess) velocity field. 
* **Pressure Correction:** This is the core of the SIMPLE algorithm. It balances the forces acting on the fluid with the mass flux across boundaries to get a more accurate pressure field, strictly fulfilling the continuity equation (mass conservation).
* **Velocity Correction:** Finally, the script corrects the initially guessed velocity field based on the newly calculated pressure gradients. 

### 3. Fighting Divergence
CFD mathematics can be highly unstable. To prevent the simulation from breaking, the code uses **relaxation coefficients** (`alpha_u`, `alpha_v`, `alpha_p`). These act as brakes on sharp gradients, slowing down the changes between iterations to ensure a smooth, stable path to convergence.

---

## Visualizations

Once the solver hits the target precision, it automatically prints the node where the flow becomes fully developed and generates a 3-part visualization figure:

1.  **Velocity Field:** A contour plot of the velocity magnitude, overlaid with solid velocity profile curves at 5 different stations to show the flow developing from a flat inlet profile to a parabola.
2.  **Pressure Field:** A contour map showing the pressure drop along the length of the pipe (relative to the outlet).
3.  **Fully Developed Comparison:** A direct comparison plot at the outlet, checking the Numerical Solver Result against the Analytical Exact Solution for fully developed laminar pipe flow.
