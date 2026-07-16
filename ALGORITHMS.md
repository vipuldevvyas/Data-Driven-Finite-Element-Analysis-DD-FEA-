# ALGORITHMS

## Algorithmic flow explicitly identified in the chat

The most compact mathematical/algorithmic sequence stated in the chat is:

`u -> B -> C -> Ke -> K`

This sequence defines the logic from interpolated displacement to assembled global stiffness system.

## Classical FEM workflow implied by the drafted Chapter 2 content

1. Define the elastic domain and boundary conditions.
2. Convert the strong form into a weak form via the principle of virtual work.
3. Discretize the domain into finite elements.
4. Approximate the displacement field with shape functions.
5. Use coordinate transformation and the Jacobian for element mapping and integration.
6. Differentiate shape functions to build the strain-displacement matrix `B`.
7. Use material law to form the constitutive matrix `C`.
8. Derive the element stiffness matrix `Ke`.
9. Assemble all `Ke` contributions into the global matrix `K` using connectivity and DOF mapping.

## Part A to Part B algorithmic break

### Part A stops at interpolation

- domain and boundary setup
- weak form
- discretization
- Jacobian
- shape functions

### Part B continues from interpolation to the system equations

- `B`
- `C`
- `Ke`
- `K`

## Validation workflow explicitly drafted for Chapter 5

For each benchmark case:

1. Define geometry and loading
2. Define material properties
3. Generate finite element mesh
4. Solve analytically if possible
5. Solve with conventional FEM
6. Solve with the developed DD-FEA model
7. Validate with Fusion 360 for 2D/3D
8. Compare displacement, strain, and stress
9. Discuss numerical accuracy and behavior

## Assembly algorithm requirement emphasized in the chat

The global assembly routine must be described as:

1. For each element, determine local DOF ordering.
2. Map local DOFs to global DOF indices using mesh connectivity.
3. Insert each entry of `Ke` into the corresponding locations of `K`.
4. Repeat for all elements until the full global matrix is assembled.

This is a preserved requirement from the chat and should not be weakened into vague prose.
