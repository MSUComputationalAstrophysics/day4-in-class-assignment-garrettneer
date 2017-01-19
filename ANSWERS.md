1. At all timesteps, the midpoint method conserves energy the best and the euler method does the worst job of conserving energy.
The Euler method seems to have linear improvement in energy change with respect to time step.
Except for the first runge-kutta point, the runge-kutta and the midpoint have 4th order improvement in energy change with respect to time step.
2. Midpoint needs about 10^2 time steps. Runge-Kutta needs about 10^3 time steps. Euler needs about 10^6 time steps
In one time step, Euler needs 4 floating point operations, midpoint needs 12 floating point operations, and runge-kutta needs 36 floating point operations.
In total, Euler needs about 4x10^6, midpoint needs 1.2x10^3, and runge-kutta needs 3.6x10^4 floating point operations.
