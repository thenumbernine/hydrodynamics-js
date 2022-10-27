[![Donate via Stripe](https://img.shields.io/badge/Donate-Stripe-green.svg)](https://buy.stripe.com/00gbJZ0OdcNs9zi288)<br>
[![Donate via Bitcoin](https://img.shields.io/badge/Donate-Bitcoin-green.svg)](bitcoin:37fsp7qQKU8XoHZGRQvVzQVP8FrEJ73cSJ)<br>

### Learning Hydrodynamics ###

![screenshot](http://christopheremoore.net/hydrodynamics/images/finite_volume_superbee_limited_sod_shock_tube_test_2d_2048x2048.png)
Screenshot is from a 2048x2048 Roe scheme with superbee flux limiter of a square Sod shock tube test.

Learning about donor cell advection, slope and flux limiters.
I currently have a Burgers and Roe solver with flux limiter implemented.

Sources: 
* Duellemond, 2009. Lecture on Hydrodynamics II http://www.mpia-hd.mpg.de/homes/dullemon/lectures/hydrodynamicsII/ 
* http://www.cfdbooks.com/cfdcodes.html 
* Toro, 1999, "Riemann Solvers and Numerical Methods for Fluid Dynamics"
* http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf



### TODO ###

1D:
- finish 'Backward Euler + Gauss-Seidel' & 'Riemann / Roe'
- more implicit options? fixed point?  some Newton/Levenberg Marquardt/BFGS stuff?
- render as a good ol fashioned white line with circles or boxes or x's or something as gridpoints
	- add a toggle for what variable is displayed
	- overlay exact solution
- HLL - incorporate flux limiter so it works at all
- AMR - isn't working unless we specifically set the flux limiter to .. (donor cell?).  the book says to use HLL, so get that working?


2D:
- custom inflow boundary conditions
	- Von Karman vortex street demo (need viscousity, i.e. Navier-Stokes, first)
- mouse input for velocity and for boundary painting
- FBO-based advection or particle simulation for streamline tracing
- implicit Gauss-Seidel ... and then add arbitrary boundaries to it

2D-GPU
- solid: get cylinder working with Roe version
- mirror boundary conditions
	- then enable hydrostatic pressure
- constant boundary options at all
- Riemann solver arbitrary boundaries
- implicit Jacobi solver

2D unstructured? 3D? 3D-GPU? 3D-unstructured? 3D-GPU-unstructured?


all:
- Entropy fix.
	- CFD codes has this in euler_solver_v0.f90, Toro describes it on p.364
- Artificial viscosity
- swappable equations of state
	- MHD
	- ADM
	- elastic solids?
- AMR

Hyperbolic Navier-Stokes:
http://deepblue.lib.umich.edu/bitstream/handle/2027.42/76470/AIAA-1991-239-186.pdf
http://aam.mathematik.uni-freiburg.de/IAM/Research/projectskr/motor/cnse.html
