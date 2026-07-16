# CODE_STRUCTURE

## Implementation context explicitly referenced in the chat

The theory is expected to match a **MATLAB DD-FEA solver**. The chat repeatedly says that the report should not become generic FEM prose detached from the implementation.

## Implementation-facing concepts referenced in the chat

- `assembleK`
- connectivity array
- DOF mapping
- insertion of local element stiffness matrices into global stiffness matrix
- progressive solver development through `1D -> 2D -> 3D`

## Implied code/report chapter alignment

The chat explicitly uses the following conceptual mapping:

- Theory progression:
  - 1D rod
  - 2D triangle / plane stress
  - 3D tetrahedron

- Implementation progression:
  - Chapter 4: 1D DD-FEA
  - Chapter 5: 2D DD-FEA / validation discussion
  - Chapter 6: 3D DD-FEA

## Report-to-code translation principle

The report should be written such that each theoretical object has an implementation counterpart:

- shape functions -> interpolation code
- Jacobian -> coordinate transformation / integration code
- `B` matrix -> strain-displacement construction
- `C` matrix -> constitutive/material routine
- `Ke` -> element stiffness routine
- global `K` -> assembly routine using connectivity and DOF indexing

## Important caution preserved from the chat

The assistant specifically warned against describing assembly as if one simply “adds all element matrices.” The documentation must reflect the actual indexing logic used in code.
