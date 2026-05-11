# Variational Discretization for Bilinear PDE-Constrained Optimal Control

A finite element implementation of a **PDE-constrained optimal control problem** with a **bilinear elliptic state equation** using **variational discretization** and **explicit segmented control projection**.

This project reproduces and studies the convergence behavior discussed in:

> A. Kröner and B. Vexler,  
> *A priori error estimates for elliptic optimal control problems with a bilinear state equation*,  
> Journal of Computational and Applied Mathematics, 230 (2009), 781–802.

---

# Problem Statement

We solve the following optimal control problem:

$$\min_{y \in H_0^1(\Omega) ,u \in L^2(\Omega} J(y,u) = \frac12 \Vert y-y_d \Vert_{L^2(\Omega)}^2 + \frac{\alpha}{2}\Vert u \Vert_{L^2(\Omega)}^2 $$

subject to the bilinear elliptic PDE

$$ -\Delta y + uy = f \quad \text{in } \quad \Omega =(0,1) $$

with homogeneous Dirichlet boundary conditions

$$ y(0)=y(1)=0 $$

and box constraints on the control

$$ 0 < a \le u(x) \le b $$

---

# Manufactured Exact Solution

The project uses a manufactured solution inspired by Example 6.2 from Kröner & Vexler:

$$ y(x)=x(1-x) $$

$$ p(x)=0.05\sin(\pi x) $$

Control obtained from the projection formula:

$$ u(x) = \mathrm{Proj}_{[a,b]} \left( \frac1\alpha y(x)p(x) \right) $$

with

$$ \alpha = 0.01, \quad a=0.1, \quad b=0.3$$

The forcing term \(f\) and desired state $y_d$ are constructed so that the exact solution satisfies the KKT system.

---

# Features

- 1D finite element discretization of state $y$
- Bilinear PDE-constrained optimization
- Variational discretization of control
- Explicit segmented control representation
- Exact polynomial integration on elements
- Exact assembly of control-dependent matrix
- Fixed-point iteration for KKT system
- L2 error computation
- Experimental Order of Convergence (EOC)
- Convergence plots
- Exact vs discrete solution visualization

---

# Numerical Method(FEM)

## State Equation

Weak formulation:

$$ \int_\Omega\nabla y \cdot \nabla v \, dx  + \int_\Omega uy\cdot v \ dx = \int_\Omega f\cdot v \ dx\qquad \qquad \forall v \in H_0^1(\Omega) $$

Discretized using linear finite elements.

## Adjoint Equation

$$ -\Delta p  + u p = (y - y_d)  \quad\text{in} \quad  \Omega $$

with boundary conditions

$$ p = 0 \quad \text{on} \quad \partial\Omega $$

Weak form:

$$ \int_\Omega (\nabla v \, \nabla p \, + \, u \, v\,p) \, dx = \int_\Omega (y - y_d) \, v \, dx \qquad \qquad \forall v \in H_0^1(\Omega) $$

## Variational Discretization

The control is **not discretized directly**.

Instead, the optimality condition is used:

$$ u(x) = \mathrm{Proj}_{[a,b]} \left( \frac1\alpha y_h(u(x)) p_h(u(x)) \right) $$

Since $y_h$ and $p_h$ are linear(P1) finite element functions:

$$ y_h p_h $$

is piecewise quadratic.

The code explicitly:
- detects active/inactive regions,
- computes roots of quadratic polynomials,
- splits elements into segments,
- integrates exactly over each segment.

---

# Convergence Results

Expected theoretical convergence:

| Variable | Expected Order |
|---|---|
| State $y_h$ | $O(h^2)$ |
| Adjoint $p_h$ | $O(h^2)$ |
| Variational Control $u_h$ | $O(h^2)$ |

The implementation computes:
- L2 errors
- EOC tables
- log-log convergence plots


# Main Components

## Exact Polynomial Integration

The implementation avoids numerical quadrature for the control matrix:

$$ \int u(x)\phi_i\phi_j \, dx $$

using:
- polynomial multiplication,
- exact symbolic integration on segments.

## Segmented Projection

Each element is split into:
- inactive quadratic region,
- lower active set,
- upper active set.

This reproduces the variational discretization idea introduced by Prof. Dr. Michael Hinze.

# Example Output

The code generates:

- EOC tables
- Convergence plots
- Exact vs discrete state
- Exact vs discrete adjoint
- Exact vs discrete control


# Mathematical Background

The optimality system consists of:

## State Equation

$$ -\Delta y^*  + u^* y^* = f \quad \text{in } \quad \Omega $$

## Adjoint Equation

$$ -\Delta p^* + u^* p^* = (y^* - y_d)  \quad\text{in} \quad  \Omega $$

## Variational Inequality

$$ (\alpha u^\ast - y^\ast p^\ast, v-u^\ast)_{L^2(\Omega)} \ge 0 $$

which leads to the projection formula

$$ u^* = \mathrm{Proj}_{[a,b]} \left( \frac{1}{\alpha} y^* (u) p^* (u) \right) $$

---

# Fixed-Point Iteration

The algorithm iteratively:

1. Assembles the control matrix
2. Solves the state equation
3. Solves the adjoint equation
4. Updates the control using projection
5. Reassembles matrices
6. Repeats until convergence

---

# Finite Element Details

The implementation uses:

- Linear finite elements (P1)
- Exact local stiffness matrices
- Exact local mass matrices
- Exact polynomial integration
- Piecewise quadratic control reconstruction

The control-dependent mass matrix is:

$$ M(u)_{ij} = \int_0^1 u(x)\phi_i(x)\phi_j(x)\,dx $$

which is assembled exactly element-by-element.

---

# Experimental Order of Convergence (EOC)

The EOC is computed using

$$\mathrm{EOC} = \frac{ \log(e_h/e_{h/2}) }{ \log(h/(h/2)) } $$

where $e_h$ denotes the numerical error.

---

# Visualizations

The project visualizes:

- Exact vs discrete state
- Exact vs discrete adjoint
- Exact vs discrete control
- Error convergence curves
- Log-log convergence rates

---

# Key Ideas Learned

This project demonstrates:

- PDE-constrained optimization
- Optimality systems
- Adjoint methods
- Variational discretization
- Finite element assembly
- Bilinear state equations
- Projection formulas
- Exact integration techniques
- Error analysis and convergence theory

---

# References

1. A. Kröner and B. Vexler,  
   *A priori error estimates for elliptic optimal control problems with a bilinear state equation*,  
   Journal of Computational and Applied Mathematics, 230 (2009), 781–802.

2. M. Hinze,  
   *Hinze, M. A Variational Discretization Concept in Control Constrained Optimization: The Linear-Quadratic Case*. Comput Optim Applic 30, 45–61 (2005). https://doi.org/10.1007/s10589-005-4559-5.

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
