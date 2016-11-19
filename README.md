# Gnedin_APC524_HW4

Solves the heat diffusion equation in C++
∂T/ ∂t = κ∇2T 
with κ = 1 on a two-dimensional domain of size 0 ≤ x ≤ π and 0 ≤ y ≤ π
With boundary conditions:
T(x, 0) = cos2 x
T(x, π) = sin2 x
T(0, y) = T(π, y) (periodic in x)
Initial condition is T=0 everywhere 

There are 3 implementations - each implementation is run on grids of size 128^2, 256^2
, and 512^2 and is run until t=0.5pi^2/kappa

serial: does not incorporate parallelism
openmp: runs with openmp using 1,2,4, and 8 threads
mpi: runs with mpi using 1,2,4,8, and 16 threads

