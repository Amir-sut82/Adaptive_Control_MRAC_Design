# Adaptive Control: MRAC Design and Robust Modifications

This repository presents a comprehensive study and implementation of **Model Reference Adaptive Control (MRAC)** techniques, including direct MRAC, indirect MRAC, integral augmentation, robust adaptive modifications, and observer-based output-feedback adaptive control.

The goal of MRAC is to adapt controller parameters online such that the uncertain plant follows the behavior of a desired stable reference model. This project investigates the performance, robustness, and limitations of different adaptive control architectures.

---

# Project Overview

The implemented tasks include:

1. **Direct MRAC** using Lyapunov-based adaptive laws
2. **Indirect MRAC** using online plant identification and certainty-equivalence control
3. **Integral augmentation** for constant disturbance rejection
4. **Robust adaptive modifications**:
   - Dead-zone modification
   - Sigma modification
   - Epsilon modification
   - Projection
5. **Output-feedback MRAC** using observer-based Loop Transfer Recovery (LTR) and adaptive augmentation
6. Comprehensive evaluation under disturbances, uncertainties, and sensor noise

The project structure follows classical adaptive control theory based on Lyapunov stability analysis. fileciteturn11file0L15-L22

---

# 1. Model Reference Adaptive Control (MRAC)

MRAC forces the plant output to follow a predefined reference model:

$$
\dot{x}_m=A_mx_m+B_mr
$$

The tracking error is:

$$
e=x-x_m
$$

The controller parameters are updated online according to the tracking error.

The main advantage of MRAC is its ability to handle matched parametric uncertainty without requiring exact knowledge of the plant parameters. However, disturbances, unmodeled dynamics, insufficient excitation, and noisy measurements can cause parameter drift. fileciteturn11file0L9-L14

---

# 2. Direct MRAC

The direct adaptive controller updates controller parameters directly using Lyapunov-based adaptation laws.

A typical control law is:

$$
u = \hat{K}_x x + \hat{K}_r r$$

The Lyapunov function is constructed as:

$$
V=e^TPe+\tilde{\theta}^T\Gamma^{-1}\tilde{\theta}
$$

The adaptive laws are designed such that:

$$
\dot{V}\leq0
$$

which guarantees boundedness of the tracking error and parameter estimates.

---

# 3. Indirect MRAC

Indirect MRAC first estimates unknown plant parameters and then designs the controller using certainty equivalence.

The identification model is:

$$
\dot{\hat{x}}=A\hat{x}+B\hat{u}+L(x-\hat{x})
$$

The estimated parameters are used to calculate the control gains.

Compared with direct MRAC, indirect MRAC can achieve improved transient performance when sufficient excitation is available.

---

# 4. Integral Augmentation

Integral augmentation introduces an additional error state to reject constant disturbances.

The integral error is:

$$
 e_I=\int e(t)dt
$$

The augmented controller improves steady-state tracking and disturbance rejection. In the project evaluation, integral augmentation significantly improves final-window RMS tracking error while maintaining similar control effort. fileciteturn11file0L215-L216

---

# 5. Robust Adaptive Modifications

Standard MRAC may suffer from parameter drift under persistent disturbances. Several robust modifications are investigated.

## Dead-zone Modification

Adaptation is stopped when the tracking error is sufficiently small:

$$
||e||<\epsilon
$$

This prevents noise-driven parameter updates.

## Sigma Modification

A leakage term is added:

$$
\dot{\hat{\theta}}=\Gamma e^TPB-\sigma\hat{\theta}
$$

which improves boundedness but introduces bias.

## Epsilon Modification

The leakage term depends on tracking error magnitude, providing stronger damping during large errors.

## Projection

Parameter estimates are constrained within a known feasible region, preventing unrealistic parameter growth.

The report discusses the trade-off between drift prevention and tracking accuracy among these methods. fileciteturn11file0L463-L467

---

# 6. Output-Feedback MRAC with LTR

When only outputs are measured, an observer is introduced:

$$
\dot{\hat{x}}=A\hat{x}+Bu+L(y-\hat{y})
$$

Loop Transfer Recovery (LTR) is used to recover desirable full-state feedback properties.

Adaptive augmentation is added to compensate matched uncertainties while keeping residual observer error and noise effects bounded. fileciteturn11file0L330-L338

---

# Simulation Evaluation

The project evaluates multiple controllers under:

- nominal tracking commands
- disturbances
- parameter variations
- sensor noise

The compared methods include direct MRAC, indirect MRAC, integral MRAC, robust modifications, and LTR-based output feedback control. fileciteturn11file378-L382

---

# Key Concepts

This repository demonstrates:

- Model Reference Adaptive Control
- Lyapunov stability analysis
- Online parameter estimation
- Robust adaptive control
- Output-feedback control
- Loop Transfer Recovery
- Disturbance rejection

---

# References

- K. S. Narendra and A. M. Annaswamy, *Stable Adaptive Systems*.
- N. T. Nguyen, *Model-Reference Adaptive Control: A Primer*.
- P. A. Ioannou and J. Sun, *Robust Adaptive Control*.
