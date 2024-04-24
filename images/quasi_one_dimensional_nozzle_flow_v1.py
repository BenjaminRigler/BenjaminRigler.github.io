# Name: Benjamin Rigler
# Year Group: LAV18
# Last edited: 14.09.2020
# This script solves a quasi one dimensional isentropic nozzle flow with the MacCormack technique.
# The example and formula is from:
# John D. Anderson, JR., "Computational Fluid Dynamics: THE BASICS WITH APPLICATIONS", McGraw-Hill
# Example 7.3
# The max CFL number is 1.1

# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
import sys

# Length of the Nozzle [m]
len_nozzle = 3

# Location of critical area [m]
location_crit_area = 1.5

# Number of cells
anzahl_zellen = 30

# Function for the area distribution [1]
area_func = lambda x: 1 + 2.2 * (x - location_crit_area) ** 2

# Function of the initial density distribution [1]
rho_func = lambda x: 1 - 0.3146 * x

# Function of the initial temperature distribution [1]
T_func = lambda x: 1 - 0.2314 * x

# Function of the initial velocity distribution [1]
v_func = lambda x, T: (0.1 + 1.09 * x) * T ** (1 / 2)

# isentropic exponent [1]
kappa = 1.4

# specific gas constant [J/kgK]
R = 287.0

# maximum number of iterations
maxit = 50000

# desired residuum
res = 10e-06

# Output
print("Start script...")
print("Print initial data...")
print("Nozzle length: ", len_nozzle, "m")
print("Location of critical area: x = ", location_crit_area, "m")
print("Number of cells: ", anzahl_zellen)
print("Isentropic Exponent: ", kappa)
print("Specific gas constant: ", R, "J/kgK")
print("Maximum number of iterations: ", maxit)
print("Desired Residuum: ", res)

# difference between nodes
delta_x = len_nozzle / anzahl_zellen

# dimensionless x coordinates
dim_x = np.linspace(0, len_nozzle, anzahl_zellen + 1)

# dimensionless area
dim_area = list(map(area_func, dim_x))

# dimensionless initial density distribution
dim_rho = list(map(rho_func, dim_x))

# dimensionless initial temperature distribution
dim_T = list(map(T_func, dim_x))

# dimensionless initial velocity distribution
dim_v = list(map(v_func, dim_x, dim_T))

# Output of initial data
print("Initial estimate of the state variables:")
print("x | Area | Density | Temperature | Velocity")
for i in range(0, anzahl_zellen + 1):
    print('%2.1f %5.3f %1.3f %1.3f %1.3f' %(dim_x[i], dim_area[i], dim_rho[i], dim_T[i], dim_v[i]))

# Input of Courant number
C = float(input("Please enter CFL Number: "))

if C <= 0:
    print("Error: Wrong input!")
    sys.exit(-1)

# Input which method for calculation of cfl number should be used
which_cfl = input("Which calculation methode of the CFL number do you want?[1 = CFL1(fast convergence) / 2 = "
                  "CFL2(slow convergence)] ")

if which_cfl == "1":
    time_calc = lambda v, T:  C * delta_x / (np.sqrt(T) + v)
elif which_cfl == "2":
    time_calc = lambda v, T: C * delta_x / (v + np.sqrt(kappa * R * T))
else:
    print("Error: Wrong Input!")
    sys.exit(-1)

# Initialising variables where all state variables along length and time are stored
T_plot = []
v_plot = []
rho_plot = []
p_plot = []

# Implement initial condition
T_plot.append(np.array(dim_T))
v_plot.append(np.array(dim_v))
rho_plot.append(np.array(dim_rho))
p_plot.append(np.multiply(dim_rho, dim_T) * R / (dim_rho[0] * dim_T[0] * R))

# Initialise new variables for calculation
# State variables
rho_new = np.array(dim_rho)
v_new = np.array(dim_v)
T_new = np.array(dim_T)

# predictor time derivatives of state variables
rho_dot = np.zeros(anzahl_zellen + 1)
v_dot = np.zeros(anzahl_zellen + 1)
T_dot = np.zeros(anzahl_zellen + 1)

# corrector time derivatives of state variables
rho_dot_cor = np.zeros(anzahl_zellen + 1)
v_dot_cor = np.zeros(anzahl_zellen + 1)
T_dot_cor = np.zeros(anzahl_zellen + 1)

# Initialise variables for residuum
res_v = []
res_rho = []
res_T = []

# Variable for iteration
it = 0

# Output
print("Start calculation...")

# Initial output of residuum
print("Residuum...")
print("Iteration | Density | Velocity | Temperature")

while True:

    # increase iteration
    it += 1

    # predictor time derivative of density
    rho_dot[1:-1] = (-np.multiply(rho_new[1:-1], (np.array(v_new[2:]) - np.array(v_new[1:-1]))) -
                     np.multiply(np.multiply(rho_new[1:-1], v_new[1:-1]), np.log(dim_area[2:]) - np.log(dim_area[1:-1]))
                     - np.multiply(v_new[1:-1], (np.array(rho_new[2:]) - np.array(rho_new[1:-1])))) / delta_x

    # predictor time derivative of velocity
    v_dot[1:-1] = (-np.multiply(v_new[1:-1], (np.array(v_new[2:]) - np.array(v_new[1:-1]))) - (np.array(T_new[2:]) -
                   np.array(T_new[1:-1]) + np.multiply(np.divide(T_new[1:-1], rho_new[1:-1]), np.array(rho_new[2:]) -
                   np.array(rho_new[1:-1]))) / kappa) / delta_x

    # predictor time derivative of temperature
    T_dot[1:-1] = (-np.multiply(v_new[1:-1], (np.array(T_new[2:]) - np.array(T_new[1:-1]))) - (kappa - 1) *
                   np.array(T_new[1:-1]) * (np.array(v_new[2:]) - np.array(v_new[1:-1]) + np.multiply(v_new[1:-1],
                   np.log(dim_area[2:]) - np.log(dim_area[1:-1])))) / delta_x

    # Calculation of time step with cfl number
    t_step = time_calc(v_new[1:-1], T_new[1:-1])

    # to have equal time steps along the geometry and stay stable, the smallest time step is taken
    t_step_min = np.min(t_step)

    # predictor state variables
    rho_temp = rho_new[1:-1] + rho_dot[1:-1] * t_step_min
    v_temp = v_new[1:-1] + v_dot[1:-1] * t_step_min
    T_temp = T_new[1:-1] + np.array(T_dot[1:-1]) * t_step_min

    # implementation of boundary condition
    v_temp = np.concatenate((2*v_temp[0] - v_temp[1], v_temp), axis=None)
    rho_temp = np.concatenate((dim_rho[0], rho_temp), axis=None)
    T_temp = np.concatenate((dim_T[0], T_temp), axis=None)

    # corrector time derivative of density
    rho_dot_cor[1:-1] = (-np.multiply(rho_temp[1:], (np.array(v_temp[1:]) - np.array(v_temp[:-1]))) -
                         np.multiply(np.multiply(rho_temp[1:], v_temp[1:]), np.log(dim_area[1:-1]) -
                         np.log(dim_area[:-2])) - np.multiply(v_temp[1:], (np.array(rho_temp[1:]) -
                         np.array(rho_temp[:-1])))) / delta_x

    # corrector time derivative of velocity
    v_dot_cor[1:-1] = (-np.multiply(v_temp[1:], (np.array(v_temp[1:]) - np.array(v_temp[:-1]))) -
                       (np.array(T_temp[1:]) - np.array(T_temp[:-1]) + np.multiply(np.divide(T_temp[1:], rho_temp[1:]),
                       np.array(rho_temp[1:]) - np.array(rho_temp[:-1]))) / kappa) / delta_x

    # corrector time derivative of temperature
    T_dot_cor[1:-1] = (-np.multiply(v_temp[1:], (np.array(T_temp[1:]) - np.array(T_temp[:-1]))) - (kappa - 1) *
                          np.multiply(T_temp[1:], (np.array(v_temp[1:]) - np.array(v_temp[:-1]) + np.multiply(v_temp[1:]
                          , np.log(dim_area[1:-1]) - np.log(dim_area[:-2]))))) / delta_x

    # corrected time derivatives
    rho_dot_average = 0.5 * (rho_dot_cor + rho_dot)
    v_dot_average = 0.5 * (v_dot_cor + v_dot)
    T_dot_average = 0.5 * (T_dot_cor + T_dot)

    # new state variables
    rho_new = rho_new + rho_dot_average * t_step_min
    v_new = v_new + v_dot_average * t_step_min
    T_new = T_new + T_dot_average * t_step_min

    # Implementation of boundary condition
    rho_new[-1:] = 2 * rho_new[-2:-1] - rho_new[-3:-2]
    v_new[-1:] = 2 * v_new[-2:-1] - v_new[-3:-2]
    v_new[0] = 2 * v_new[1] - v_new[2]
    T_new[-1:] = 2 * T_new[-2:-1] - T_new[-3:-2]

    # Calculation of pressure with state equation for ideal gases
    p_plot.append(R * T_new * rho_new / (dim_rho[0] * dim_T[0] * R))

    # Store the new state variables
    rho_plot.append(rho_new)
    v_plot.append(v_new)
    T_plot.append(T_new)

    # RMS values of corrected time derivatives (is used as Residuum)
    rho_dot_rms = np.sqrt(np.power(rho_dot_average, 2).sum(axis=None))
    v_dot_rms = np.sqrt(np.power(v_dot_average, 2).sum(axis=None))
    T_dot_rms = np.sqrt(np.power(T_dot_average, 2).sum(axis=None))

    # store the Residuum to plot further on
    res_rho.append(rho_dot_rms)
    res_v.append(v_dot_rms)
    res_T.append(T_dot_rms)

    # Output of the residuum
    print('%2.1i %s %1.7f %s %1.7f %s %1.7f' % (it, " | ", rho_dot_rms, " | ", v_dot_rms,
                                                      " | ", T_dot_rms))

    # Stop condition
    if rho_dot_rms <= res and v_dot_rms <= res and T_dot_rms <= res:
        print("Calculation successful...")
        break
    if it == maxit:
        print("Calculation stopped! Max Iterations reached...")
        break

# Output of the different plot options
print("Choose your plot...")
print("g = Geometry")
print("t = State variables as a function of time")
print("l = State variables as a function of length")
print("r = Residuum")
print("e = End script")

# Input loop
while True:

    # Input
    which_plot = input("...: ")

    # State variables as function of the time
    if which_plot == "t":

        # Input of the x-coordinate
        x = float(input("Please enter x-coordinate: "))

        # initialise variables to store interpolated data
        rho_interp = []
        v_interp = []
        T_interp = []
        p_interp = []

        # interpolate the state variables in the chosen x-coordinate
        for i in rho_plot:
            rho_interp.append(np.interp(x, dim_x, i))

        for i in v_plot:
            v_interp.append(np.interp(x, dim_x, i))

        for i in T_plot:
            T_interp.append(np.interp(x, dim_x, i))

        for i in p_plot:
            p_interp.append(np.interp(x, dim_x, i))

        # plot
        x_i = np.linspace(0, it, it + 1)

        fig1, (ax1_1, ax1_2, ax1_3, ax1_4) = plt.subplots(4, sharex=True)
        fig1.suptitle("State variables as function of time")
        ax1_1.plot(x_i, rho_interp, 'r')
        ax1_2.plot(x_i, v_interp, 'b')
        ax1_3.plot(x_i, T_interp, 'g')
        ax1_4.plot(x_i, p_interp, 'y')
        ax1_1.xaxis.grid(color='gray')
        ax1_1.yaxis.grid(color='gray')
        ax1_2.xaxis.grid(color='gray')
        ax1_2.yaxis.grid(color='gray')
        ax1_3.xaxis.grid(color='gray')
        ax1_3.yaxis.grid(color='gray')
        ax1_4.xaxis.grid(color='gray')
        ax1_4.yaxis.grid(color='gray')
        ax1_1.set(ylabel="$\\rho$/$\\rho_0$")
        ax1_2.set(ylabel="M")
        ax1_3.set(ylabel="$T/T_0$")
        ax1_4.set(xlabel="i", ylabel="$p/p_0$")
        plt.show()

    # state variables as function of length
    elif which_plot == "l":

        fig2, (ax2_1, ax2_2, ax2_3, ax2_4) = plt.subplots(4, sharex=True)
        fig2.suptitle("State variables as function of length")
        ax2_1.plot(dim_x, rho_plot[-1], 'r')
        ax2_2.plot(dim_x, v_plot[-1], 'b')
        ax2_3.plot(dim_x, T_plot[-1], 'g')
        ax2_4.plot(dim_x, p_plot[-1], 'y')
        ax2_1.xaxis.grid(color='gray')
        ax2_1.yaxis.grid(color='gray')
        ax2_2.xaxis.grid(color='gray')
        ax2_2.yaxis.grid(color='gray')
        ax2_3.xaxis.grid(color='gray')
        ax2_3.yaxis.grid(color='gray')
        ax2_4.xaxis.grid(color='gray')
        ax2_4.yaxis.grid(color='gray')
        ax2_1.set(ylabel="$\\rho$/$\\rho_0$")
        ax2_2.set(ylabel="M")
        ax2_3.set(ylabel=r"$T/T_0$")
        ax2_4.set(xlabel="x/L", ylabel=r"$p/p_0$")
        plt.show()

    # residuum
    elif which_plot == 'r':

        x_i = np.linspace(0, it-1, it)

        fig3 = plt.figure()
        plt.plot(x_i, res_v, 'r', label="Residuum Velocity")
        plt.plot(x_i, res_rho, 'b', label="Residuum Density")
        plt.plot(x_i, res_T, 'g', label="Residuum Temperature")
        plt.title("Residuum")
        plt.xlabel("Iterations")
        plt.ylabel("RMS")
        plt.legend(loc='upper right')
        plt.grid()
        plt.show()

    # end script
    elif which_plot == 'e':
        break

    # geometry
    elif which_plot == "g":

        area = list(map(area_func, dim_x))

        fig4 = plt.figure()
        plt.plot(dim_x, area, 'gray')

        area = np.multiply(area, -1)

        plt.plot(dim_x, area, 'gray')
        plt.title("Geometry")
        plt.xlabel("x/L")
        plt.ylabel(r"$A/A_{critical}$")
        plt.grid()
        plt.show()

    else:
        print("Error: Wrong input!")


