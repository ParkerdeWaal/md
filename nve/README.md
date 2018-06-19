Initial NVE practice system with the following potential energy function

```fortran
V(x,y) = 4(x^2+y^2-1)^2y^2 - e^{-4((x-1)^2+y^2)} - e^{-4((x+1)^2+y^2)} + e^{8(x-1.5)} + e^{8(x-+.5)} + e^{-4(y+0.25)} + 0.2e^{-8x^2}
```
To compile and run:
`gfortran -O3 nve.f90`
`./a.out > output`

Output file contains time, hamiltonian, x, and y coordinates.

to visualize output in gnuplot:
`load'surf.gnu`
`replot 'output' u 3:4:2 w p`

