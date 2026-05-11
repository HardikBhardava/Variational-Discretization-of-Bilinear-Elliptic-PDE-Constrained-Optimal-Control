# PDE-Constrained Optimal Control:
# Constant (P0) vs Variational Control Discretization

This project implements and compares two finite element approaches for a **bilinear elliptic PDE-constrained optimal control problem**:

1. **Piecewise constant (P0) control discretization**
2. **Variational discretization with segmented control reconstruction**

The implementation studies:
- finite element discretization,
- optimality systems,
- control projection formulas,
- convergence behavior,
- experimental orders of convergence (EOC),
- exact vs discrete solution comparison.

The project is inspired by the variational discretization framework introduced by Prof.Dr. Michael Hinze and the error analysis of Kröner & Vexler.

---

# PDE-Constrained Optimization Problem

We solve the optimal control problem

$$\min_{y \in H_0^1(\Omega) ,u \in L^2(\Omega} J(y,u) = \frac12 \Vert y-y_d \Vert_{L^2(\Omega)}^2 + \frac{\alpha}{2}\Vert u \Vert_{L^2(\Omega)}^2 $$

subject to the bilinear elliptic PDE

$$ -\Delta y + uy = f \quad \text{in } \quad \Omega =(0,1) $$

with homogeneous Dirichlet boundary conditions

$$
y(0)=y(1)=0
$$

and box constraints on the control

$$ 0 < a \le u(x) \le b $$

---

# Optimality System

The first-order optimality system consists of:

## State Equation

$$ -\Delta y + uy = f \quad \text{in } \quad \Omega =(0,1) $$

## Adjoint Equation

$$ -\Delta p + up = y - y_d \quad \text{in } \quad \Omega =(0,1) $$

## Variational Inequality

$$ (\alpha u^\ast - y^\ast p^\ast, v-u^\ast)_{L^2(\Omega)} \ge 0 $$

which leads to the projection formula

$$ u^* = \mathrm{Proj}_{[a,b]} \left( \frac{1}{\alpha} y^* (u) p^* (u) \right) $$

---

# Manufactured Exact Solution

The code uses a manufactured exact solution:

## Exact State

$$
y(x)=x(1-x)
$$

## Exact Adjoint

$$
p(x)=0.05\sin(\pi x)
$$

## Exact Control

$$
u(x)
=
\mathrm{Proj}_{[a,b]}
\left(
\frac{1}{\alpha} y(x)p(x)
\right)
$$

with

$$
\alpha=0.01,
\quad
a=0.1,
\quad
b=0.3
$$

The forcing term $f$ and desired state $ y_d $ are constructed so that the exact solution satisfies the full KKT system.

---

# Numerical Methods(FEM)

The project compares two different control discretization strategies.

---

# 1. Piecewise Constant (P0) Control

The control is discretized elementwise:

$$
u_h|_K = \text{constant}
$$

The algorithm:
1. assembles FEM matrices,
2. solves state equation,
3. solves adjoint equation,
4. updates control elementwise,
5. applies projection onto \([a,b]\),
6. repeats until convergence.

The control update is based on the averaged quantity

$$
u_h
=
\mathrm{Proj}_{[a,b]}
\left(
\frac{1}{\alpha h}
\int_K y_hp_h\,dx
\right)
$$

Expected convergence:

| Variable | Expected Order |
|---|---|
| State | \(O(h^2)\) |
| Adjoint | \(O(h^2)\) |
| Control | \(O(h)\) |

---

# 2. Variational Discretization

In the variational approach:
- the control is **not directly discretized**,
- only state and adjoint are discretized using P1 FEM.

The control is reconstructed from

$$
u_h(x)
=
\mathrm{Proj}_{[a,b]}
\left(
\frac{1}{\alpha} y_h(x)p_h(x)
\right)
$$

Since $y_h$ and $p_h$ are linear finite element functions,

$$
y_h(x)p_h(x)
$$

is piecewise quadratic.

The implementation:
- reconstructs quadratic polynomials,
- computes roots analytically,
- detects active/inactive regions,
- splits elements into segments,
- performs exact polynomial integration.

Expected convergence:

| Variable | Expected Order |
|---|---|
| State | \(O(h^2)\) |
| Adjoint | \(O(h^2)\) |
| Control | \(O(h^2)\) |

---

# Main Features

- 1D finite element method (P1 FEM)
- Bilinear PDE-constrained optimization
- Fixed-point iteration
- Variational discretization
- Piecewise constant control discretization
- Exact polynomial integration
- Exact segmented control assembly
- Active/inactive set reconstruction
- L2 error computation
- Experimental Order of Convergence (EOC)
- Convergence plots
- Exact vs numerical solution comparison

---

# Finite Element Formulation

The weak formulation of the state equation is:

$$ \int_\Omega\nabla y \cdot \nabla v \ dx  + \int_\Omega uy\cdot v \ dx = \int_\Omega f\cdot v \ dx\qquad \qquad \forall v \in H_0^1(\Omega) $$

using linear finite element basis functions.

---

# Control-Dependent Matrix

The bilinear term produces the matrix

$$
C(u)_{ij}
=
\int_0^1 u(x)\phi_i(x)\phi_j(x)\,dx
$$

which is assembled:
- exactly for the variational method,
- elementwise constant for the P0 method.

---

# Exact Polynomial Integration

The variational discretization uses:
- polynomial multiplication,
- exact integration over subsegments,
- exact local matrix assembly.

This avoids quadrature errors in the control reconstruction.

---

# Segmented Control Reconstruction

Each finite element is partitioned into:
- lower active region,
- inactive quadratic region,
- upper active region.

The code analytically computes:
- roots of quadratic polynomials,
- transition points between active sets,
- exact integrals on each segment.

---

# Experimental Order of Convergence (EOC)

The EOC is computed using

$$
\mathrm{EOC}
=
\frac{
\log(e_h/e_{h/2})
}{
\log(h/(h/2))
}
$$

where $e_h$ is the numerical error.

---

# Output

The code generates:

- EOC tables
- State convergence plots
- Adjoint convergence plots
- Control convergence plots
- Exact vs discrete state
- Exact vs discrete adjoint
- Exact vs discrete control
- 
---

# Numerical Observations

The implementation confirms the theoretical behavior:

| Method | Control Convergence |
|---|---|
| P0 constant control | \(O(h)\) |
| Variational discretization | \(O(h^2)\) |

while state and adjoint achieve approximately

$$
O(h^2)
$$

for both methods.

---

# References

1. Michael Hinze,  
   *Hinze, M. A Variational Discretization Concept in Control Constrained Optimization: The Linear-Quadratic Case*. Comput Optim Applic 30, 45–61 (2005). https://doi.org/10.1007/s10589-005-4559-5.

2. A. Kröner and B. Vexler,  
   *A priori error estimates for elliptic optimal control problems with a bilinear state equation*,  
   Journal of Computational and Applied Mathematics, 230 (2009), 781–802. :contentReference[oaicite:0]{index=0}

3. Fredi Tröltzsch,  
   *Optimal Control of Partial Differential Equations*, Springer.

---

# Author

Hardik Bhardava

Research interests:
- PDE-constrained optimization
- Finite element methods
- Scientific machine learning
- Physics-informed neural networks
- Numerical analysis
