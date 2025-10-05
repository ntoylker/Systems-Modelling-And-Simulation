# 🖥️ Systems Modelling and Simulation – University Project Repository

This repository contains my work for the **Systems Modelling and Simulation** university course, spring 2025.  
It includes solved projects, simulations, and scripts for system modeling and parameter estimation, written in **MATLAB**.

---

## Project 1: Estimation of Unknown Parameters – Least Squares Method  
**Date:** 20 March 2025  

**System:** A simple pendulum with input torque.  
Linearized equation:  
`
m L^2 \ddot{q}(t) + c \dot{q}(t) + m g L q(t) = u(t)
`
where `q(t)` is the angular displacement, `m` mass, `L` length, `c` damping coefficient, `g` gravity, `u(t)` control input.

**Tasks & Methods:**
- **State-space representation & transfer function**  
  - Derived `ẋ(t) = A x(t) + B u(t)`  
  - Output is `q(t)`  

- **System simulation**  
  - Solved ODE using MATLAB with zero initial conditions  
  - Input: `u(t) = A0 sin(ω t)`  
  - Parameters: `m=0.75`, `L=1.25`, `c=0.15`, `g=9.81`, `A0=4`, `ω=2`  
  - Simulated for 20 sec with high precision (Δt < 10⁻³ sec)  
  - Generated plots of system states  

- **Parameter estimation using Least Squares**  
  - Using measured `x(t)` and `u(t)` to estimate `m`, `L`, `c`  
  - Compared system response with estimated parameters vs. original  
  - Plotted `q(t)`, estimated `q̂(t)` and error `e_q(t)`  

- **Partial measurements**  
  - Applied Least Squares using only `q(t)` and `u(t)`  
  - Analyzed results  

- **Robustness & sensitivity analysis**  
  - Added Gaussian noise to measurements and recomputed parameter estimates  
  - Studied effect of sampling period `Ts` on estimation accuracy  
  - Studied effect of input amplitude `A0` on estimation accuracy  

**Methods Used:** MATLAB simulations, ODE solvers, Least Squares estimation, noise analysis, sensitivity studies.

---

## Project 3: Real-Time Robust Methods – Projection Methods  
**Date:** 23 May 2025  

**System:** Linear system  
`
ẋ(t) = A x(t) + B u(t)
`
where `x(t)` ∈ ℝ² is the state, `u(t)` ∈ ℝ is the input, `A` and `B` unknown but bounded matrices.

**Tasks & Methods:**
- **Real-time parameter estimation algorithm**  
  - Designed algorithm to estimate unknown matrices `A` and `B`  
  - Input `u(t)` chosen to excite the system  
  - Evaluated stability of estimation  
  - Plotted `x(t)`, estimated `x̂(t)` and error `e_x(t)`  
  - Plotted estimates `Â(t)` and `B̂(t)`  

- **Modeling and bias handling**  
  - Introduced a bounded bias `ω(t)` in the system: `||ω(t)|| ≤ ω̄`  
  - Developed modeling approach for bias  
  - Applied real-time estimation algorithm under biased input  
  - Studied effect of bias magnitude `ω̄` on estimation accuracy  

**Methods Used:** Real-time estimation, projection methods, stability analysis, MATLAB simulations, sensitivity analysis.

---

**Technologies and Tools Used:**  
- MATLAB & Simulink for simulations  
- Linear system modeling, Least Squares estimation, real-time projection methods
