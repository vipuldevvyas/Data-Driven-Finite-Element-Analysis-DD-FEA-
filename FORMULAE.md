# FORMULAE

## Core formula chain preserved from the chat

- Displacement interpolation:
  - `u = N d_e`
- Strain-displacement relation:
  - `epsilon = B d_e`
- Constitutive relation:
  - `sigma = C epsilon`
- Combined relation:
  - `sigma = C B d_e`
- Element stiffness concept:
  - `Ke = integral(B^T C B dOmega)`
- Global assembly concept:
  - local `Ke` matrices are inserted into global `K` through connectivity and DOF mapping

## Governing equations referenced in drafted Chapter 2 Part A

### Domain and boundaries

- body domain: `Omega`
- full boundary: `Gamma`
- displacement boundary: `Gamma_u`
- traction boundary: `Gamma_t`

### Linear elasticity equilibrium

- `div(sigma) + b = 0`

### Boundary conditions

- displacement condition:
  - `u = u_bar` on `Gamma_u`
- traction condition:
  - `sigma n = t_bar` on `Gamma_t`

## Weak form / virtual work relations

The chat’s draft explicitly moves from strong form to weak form using virtual displacement and the divergence theorem. The resulting key idea is:

- internal virtual work = external virtual work

This is presented as the mathematical basis of FEM.

## Shape-function interpolation relations referenced in the chat

### General interpolation

- `u = N d`

### 1D bar element shape functions

- `N1 = 1 - x/L`
- `N2 = x/L`

### 2D triangle in natural coordinates

- `N1 = 1 - xi - eta`
- `N2 = xi`
- `N3 = eta`

### 3D tetrahedron in natural coordinates

- `N1 = 1 - xi - eta - zeta`
- `N2 = xi`
- `N3 = eta`
- `N4 = zeta`

## Jacobian-related relations

The draft explicitly identifies:

- physical-to-natural coordinate transformation
- Jacobian matrix `J`
- determinant condition:
  - `|J| > 0`

## B-matrix content preserved from the chat

### 1D bar

- strain:
  - `epsilon = du/dx`
- B matrix:
  - `[-1/L, 1/L]`

### 2D triangular element

The chat explicitly says the full `3 x 6` matrix should be written using:

- `dNi/dx`
- `dNi/dy`

and strain components:

- `epsilon_x`
- `epsilon_y`
- `gamma_xy`

### 3D tetrahedral element

The chat explicitly says the full `6 x 12` tetrahedral `B` matrix is written in complete form using derivatives of:

- `N1`
- `N2`
- `N3`
- `N4`

with respect to:

- `x`
- `y`
- `z`

## C-matrix content preserved from the chat

### 1D constitutive relation

Not fully written in visible excerpt, but explicitly recommended as part of the final `1D -> 2D -> 3D` progression for Section `2.8`.

### 2D plane stress matrix

The chat explicitly requires that the plane stress constitutive matrix be written completely.

### 3D isotropic elasticity matrix

The chat explicitly requires that the full `6 x 6` constitutive matrix be written completely.

### Lamé-parameter form

The visible draft includes:

- first Lamé parameter `lambda`
- shear modulus `mu`

and states that the constitutive matrix can also be written in Lamé form.

## Most important formula-related style rule

Every major equation should be followed by `1-2` explanatory lines in the final report.
