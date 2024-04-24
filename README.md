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
![](images/advdiff/diff_meshComp.png) |
![](images/advdiff/diff_lineplot.png) |

Convection             |
:-------------------------:|
![](images/advdiff/swirl.gif) |

# Navier-Stokes using Finite-Volumes and the Simple-algorithm
This project solves the 2D unsteady and steady Navier_Stokes equations on a staggered cartesian grid. The convection term can be discretized using a upwind, central or hybrid differencing scheme. The diffusion term is discretized using central differences and the time derivative is discretized using the Adams_Moulton scheme with variable time stepping. The results are stored in the .vtk file format, which makes post processing in ParaView possible. Next a comparison of the velocity with data from Ghia et al. is shown for a steady state solution and a subsequent animation of u-velocity component with streamlines for the unsteady solution. All simulations have Re=100.

Steady-state            |
:-------------------------:|
![](images/simple/post_Re100_u_v2.png) |
![](images/simple/post_Re100_v_v2.png) |


Steady-state            |
:-------------------------:|
![](images/simple/animation_re100.gif) |

# GPU-Accelerated Euler-Solver
This project is based on another project done during the course Numerical Methods for Fluid Mechanics at TU Wien. The original goal was to implement the Jameson-algorithm to solve the Euler equations over canal with a circular bump in Python. However this code didn't fully work, so I decided to give it another run in my free time. I reimplemented the solver in C++ and could fix the issue. Furthermore I tried to improve performance by using GPU computing (CUDA). It seems that the code running on the GPU has a higher floating point error than the code running on the CPU, which can be fixed by increasing artificial diffusion. Next some results are shown, however the project is not fully finished yet, so proper validation and verification with literature and a performance benchmark is still pending.

Flow over circular bump           |
:-------------------------:|
![](images/euler/bump.gif) |

Flow in nozzle           |
:-------------------------:|
![](images/euler/nozzle.gif) |
