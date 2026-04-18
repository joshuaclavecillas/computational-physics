# Notebook Summaries

## 1D_molecular_dynamics_velocity_verlet.ipynb
This notebook implements a 1D molecular dynamics simulation using the Velocity Verlet integration algorithm to model particles interacting via a Lennard-Jones potential in a periodic box. The simulation's energy conservation is tested, and the effects of timestep size and temperature-based initial velocities on the system's behavior are investigated.

## 2D Rock-Paper-Scissors Territorial Walkers.ipynb
This notebook simulates a 2D territorial battle between three species: Rock, Paper, and Scissors, where walkers of each species move randomly on a grid. When walkers of different species meet, they convert territory based on RPS rules, with added stochastic elements to explore how initial conditions and randomness influence the long-term dominance or coexistence of the species.

## 2D_molecular_dynamics_velocity_verlet.ipynb
This notebook simulates the two-dimensional motion of particles in a box with reflecting walls, using the Velocity Verlet method for time integration. The particles interact through a repulsive Lennard-Jones potential, and the simulation tracks energy conservation and the system's temperature over time.

## 3D Random Walk Simulation.ipynb
This notebook simulates and analyzes three-dimensional random walks where each step is normalized to a unit length. It visualizes multiple independent walks, calculates statistical properties like mean displacement, and validates the randomness by checking for correlations between steps.

## Earth_Sun_Orbit_RK2.ipynb
This notebook simulates the Earth's orbit around the Sun in two dimensions using the second-order Runge-Kutta (midpoint) method to integrate the equations of motion. It uses realistic orbital parameters for Earth, including its eccentricity, and analyzes the conservation of energy and angular momentum to validate the simulation's accuracy.

## Example 2.7 (computational).ipynb
This notebook contains two small Python code snippets. The first calculates and prints a sequence of numbers based on a recursive formula, while the second converts a pair of Cartesian coordinates (x, y) into polar coordinates (r, theta).

## Hw (2.4, 2.5, 2.9).ipynb
This notebook contains solutions to three homework problems. The first problem calculates the time dilation for a spaceship traveling at 99% the speed of light, the second computes the transmission and reflection probabilities of an electron encountering a potential step, and the third attempts to calculate the Madelung constant for a sodium chloride crystal.

## Metropolis_algorithm_demo.ipynb
This notebook demonstrates the Metropolis algorithm by simulating a 1D Ising model of spin configurations. It implements the spin-flip update rule, where spins are randomly flipped and the new configuration is accepted or rejected based on the change in energy and the system's temperature. The simulation tracks the evolution of magnetization and energy, and visualizes the spin configurations over time to show how the system reaches equilibrium.

## Pair work 2.ipynb
This notebook contains solutions for two exercises. The first exercise performs and verifies simple array arithmetic using NumPy, while the second exercise calculates the binding energy of an atomic nucleus using the semi-empirical mass formula and explores how the binding energy per nucleon varies with the mass number and atomic number.

## Pentagon_gasket_starter.ipynb
This notebook implements the "chaos game" to generate a fractal pattern based on a regular pentagon. It starts with a random point inside the pentagon and iteratively moves it halfway toward a randomly chosen vertex, plotting the resulting points to reveal a complex, self-similar structure. The notebook also includes an animation that zooms in on one of the vertices to show the fractal's intricate details.

## Point_Charge_Spring_Euler_Simulation.ipynb
This notebook simulates the motion of a point charge attached to a spring using the explicit Euler method. It starts with a simple harmonic oscillator, analyzes its energy conservation, and then introduces additional forces like random wind gusts and spatial drag to observe their effects on the system's dynamics and energy. The notebook also explores anharmonic oscillations by modifying the potential to a power-law form and investigates the numerical stability of the Euler method for these nonlinear systems.

## Projectile_Motion_with_Drag_Simulation.ipynb
This notebook simulates the trajectory of a projectile, such as a ball, considering both gravity and air resistance (drag). It uses the Euler method to numerically solve the equations of motion and visualizes the resulting trajectory, velocity, and speed over time. The simulation also compares the numerical results with the analytical solution for projectile motion without drag to quantify the effects of air resistance on the flight path.

## Projectile_with_Drag.ipynb
This notebook simulates projectile motion with linear air drag using `scipy.integrate.solve_ivp` to solve the differential equations of motion. It compares the trajectory with and without drag, validates the numerical solution against the analytical solution, and explores how the trajectory changes with different drag coefficients and launch angles. The notebook also includes an analysis to find the optimal launch angle for maximum range under the influence of drag.

## RPS_circular_territories.ipynb
This notebook simulates a territorial battle between Rock, Paper, and Scissors on a grid with a different initial setup: the three species start in circular sectors, creating a Y-shaped symmetry. The simulation uses the same walker mechanics as previous versions, where walkers move randomly and convert territory based on RPS rules. The notebook visualizes the evolution of the territories and includes an animation to show how the initial circular boundaries break down and evolve into complex, chaotic patterns.

## RPS_random_spawn.ipynb
This notebook simulates a territorial battle between Rock, Paper, and Scissors where the walkers of each species are randomly spawned across the grid at the beginning of the simulation. It uses the same core mechanics as the other RPS simulations, with walkers moving and converting territory based on the RPS rules. The notebook explores how the initial random distribution of walkers influences the evolution of the system, including the emergence of dominant species and the formation of complex spatial patterns.

## RPS_vertical_thirds.ipynb
This notebook simulates a territorial battle between Rock, Paper, and Scissors on a grid where the initial territories are divided into three vertical sections. The simulation uses the same walker mechanics as previous versions, with walkers moving randomly and converting territory based on RPS rules. The notebook visualizes the evolution of the territories and includes an animation to show how the initially straight boundaries become chaotic and complex over time.

## Random Walk Simulation.ipynb
This notebook simulates a 2D random walk with normalized steps, visualizing the paths and calculating displacement statistics over many trials. It also verifies the core assumption of a random walk by showing that cross-step correlations average to zero for a large number of steps, and confirms the theoretical relationship where the root-mean-square displacement scales with the square root of the number of steps.

## Sierpinski_gasket_starter.ipynb
This notebook demonstrates the "chaos game" method for generating the Sierpiński gasket fractal. It iteratively plots points by starting with an arbitrary point inside a triangle and repeatedly moving it halfway toward a randomly selected vertex, revealing the classic fractal structure. The notebook also creates an animation that zooms into one of the vertices, highlighting the self-similar nature of the gasket at different scales.

## Uranus_Orbit_Sun_RK45.ipynb
This notebook simulates the orbit of Uranus around the Sun using the adaptive RK45 numerical integrator. It first models a simple two-body system to establish a baseline orbit and then introduces a Neptune-like body to demonstrate how mutual gravitational interactions perturb the orbit, calculating the resulting change in Uranus's orbital period.

## box_counting_barca_crest.ipynb
This notebook applies the box-counting algorithm to an image of the FC Barcelona crest to estimate its fractal dimension. It analyzes both a foreground mask and an edge map of the crest, calculating the dimension from the slope of a log-log plot of box counts versus box size. The notebook also visualizes the box-fitting process, showing how the image is covered by grids of different sizes to measure its complexity.

## calculate_4x2.ipynb
This is a very simple notebook that contains a single calculation. It computes the product of 4 and 2.

## capacitor_problem.ipynb
This notebook uses the Jacobi relaxation method to numerically solve for the electrostatic potential and electric field in various configurations. It begins with a 2D parallel-plate capacitor, progresses to a 3D charged torus, and culminates in a system containing both a charged torus and an oppositely charged sphere, complete with an animation of the resulting electric field streamlines.

## capacitor_problem_with_pseudocode.ipynb
This notebook provides a detailed, step-by-step guide to solving Laplace's equation for a 2D capacitor, starting with pseudocode for the Jacobi relaxation method. It then systematically compares the performance of the Jacobi method against the Gauss-Seidel and Successive Over-Relaxation (SOR) methods, demonstrating how more advanced iterative solvers can significantly improve convergence speed.

## exercise 2.12, 3.1, 3.3.ipynb
This notebook contains solutions to three distinct programming exercises. It includes a prime number generator, a data analysis script that plots and smooths sunspot activity data using a running average, and a visualization tool that loads a 2D dataset from "stm.txt" and displays it as a color-mapped image.

## homework 2(3.7 and 4.4).ipynb
This notebook presents solutions to two homework problems: generating the Mandelbrot set and performing numerical integration. The first part visualizes the Mandelbrot set by iterating a complex function, while the second part calculates the integral of a semicircle using both a custom trapezoidal rule implementation and NumPy's built-in function for comparison.

## maxwell_1D_box_fourier_sol.ipynb
This Jupyter notebook file is currently empty and does not contain any code or content.

## nested_boxes_field_calculation.ipynb
This notebook calculates the electric potential and field for a nested box capacitor using the successive over-relaxation (SOR) method. It visualizes the results by plotting the potential as a colormap with equipotential lines and the electric field magnitude with overlaid streamlines, clearly showing the field concentration around the corners of the inner box.

## pair work 3.ipynb
This notebook tackles several numerical methods exercises, starting with a factorial calculation that highlights data type limits and an exploration of numerical differentiation errors using the forward difference formula. It then implements a more accurate central difference method to compute a derivative and concludes by applying numerical gradients to an altitude map to generate a relief shading visualization of the terrain.

## radioactive_decay_simulation.ipynb
This notebook simulates radioactive decay using a Poisson distribution, first with a fixed decay constant and then with a randomly varying one to compare the outcomes. It visualizes the exponential decay on a log scale and analyzes the relationship between the decay constant, initial particle count, and decay rate, confirming the expected linear relationship between the logarithm of the particle count and the logarithm of the change in particle count.

## simple_pendulum_rk4_phase_diagram.ipynb
This notebook simulates a pendulum-like system using a 4th-order Runge-Kutta (RK4) integrator to solve the governing differential equation. It first visualizes the time-series and phase diagram for a single trajectory, then generates a comprehensive phase portrait by overlaying the trajectories from a wide range of initial conditions.

## stable_population_logistic_map.ipynb
This notebook models population dynamics using the logistic map, analyzing its fixed points and their stability by numerically estimating the function's derivative. It visualizes the population's evolution over many generations for different growth rates (μ), clearly demonstrating the transition from stable points to period-doubling and chaos.

## test_corrections.ipynb
This notebook attempts to simulate projectile motion under gravity using the Euler method. The code initializes parameters for the simulation but contains significant syntax and logical errors in its main loop, preventing it from running correctly.

## week 5 work.ipynb
This notebook explores the properties of a Linear Congruential Generator (LCG), demonstrating its periodic nature and the non-random correlations that appear when plotting successive pairs of generated numbers. For contrast, it also visualizes pairs from NumPy's random number generator to highlight the difference between a simple LCG and a more sophisticated pseudo-random number algorithm.
