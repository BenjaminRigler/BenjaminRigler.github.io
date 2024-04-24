# Quasi-one-dimensional nozzle flow
The goal of this project is to solve instationary quasi-one-dimensional flow in a nozzle using the MacComack-method. The method is described in John D. Anderson, JR., "Computational Fluid Dynamics: THE BASICS WITH APPLICATIONS", McGraw-Hill

![](images/nozzle_flow/geometry.png)
![](images/nozzle_flow/results.png)

# Navier-Stokes using the artificial compressibility
This project solves the Navier-Stokes equations using the artificial compressibility method. A simple pre- and post-processor is also implemented within the program. The results for the single lid driven cavity problem are shown next. Upon comparison with results from Ghia et al. the results seem correct.

u-velocity             |
:-------------------------:|
![](images/arti_comp/u_contour.png)  |
![](images/arti_comp/u_line.png)  |

v-velocity             |
:-------------------------:|
![](images/arti_comp/v_contour.png) |
![](images/arti_comp/v_line.png) |

# Unsteady convection-diffusion equation on non-orthogonal
The goal of this project was to experience the challenges of dealing with non-orthogonal meshes. For that a solver for the 2d unsteady convection-diffusion equation was implemented using the Finite-Volume method. For discretisation of the convection term HO schemes are supported. The diffusion term is discretizes using a central difference scheme and the time derivative is discretized using the implizit Euler method. To keep the generation of the mesh as flexible as possible the format implemented by OpenFOAM is used. Next the solution of a diffusion problem is compared between a cartesian and an unstructured triangular mesh. It can be observed that the code can deal with non-orthogonal problems. Furthermore a convection problem was tested. Here the solver does not converge for the unstructured mesh, which is a significant issue of the code to be fixed.

Diffusion             |
:-------------------------:|
![](images/advconv/diff_meshComp.png)
![](images/advconv/diff_lineplot.png)

Convection             |
:-------------------------:|
![](images/advconv/swirl.gif)

