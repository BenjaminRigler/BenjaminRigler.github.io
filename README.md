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
