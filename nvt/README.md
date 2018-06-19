Canonical ensemble (NVT) for a single particle system using a symplectic integrator for Langevin's equation. The the following potential energy function is used here:

```fortran
V(x,y) = 4(x^2+y^2-1)^2y^2 - e^{-4((x-1)^2+y^2)} - e^{-4((x+1)^2+y^2)} + e^{8(x-1.5)} + e^{8(x-+.5)} + e^{-4(y+0.25)} + 0.2e^{-8x^2}
```

To compile and run:
`gfortran -O3 nvt.f90 random_numbers.f`
`./a.out > output`

Output file contains time, hamiltonian, kinetic energy, x, and y coordinates.

to visualize output in gnuplot:
`load'surf.gnu`
`replot'output' u 4:5:2 w p`

