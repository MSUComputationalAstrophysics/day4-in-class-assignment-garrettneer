1. At all timesteps, the runge-kutta method conserves energy the best and the euler method does the worst job of conserving energy.
The Euler method seems to have linear improvement in energy change with respect to time step.
The midpoint seems to have 3rd order improvement in energy change with respect to time step.
The runge-kutta seems to have 4th order improvement in energy change with respect to time step.
2. Runge-Kutta needs about 50 time steps. Midpont needs about 200 time steps. Euler needs about 10^6 time steps
In one time step, Euler needs 4 floating point operations, midpoint needs 12 floating point operations, and runge-kutta needs 36 floating point operations.
In total, Euler needs about 4x10^6, midpoint needs 2400, and runge-kutta needs 1800 floating point operations.
Even though runge-kutta has better improvement with respect to time step compared to midpoint, it does not save much time.
