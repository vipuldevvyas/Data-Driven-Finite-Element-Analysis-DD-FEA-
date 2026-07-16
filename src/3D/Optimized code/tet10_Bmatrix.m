function [B_all, V_all, GP_all] = tet10_Bmatrix(coords)

% Tet10 B-matrix computation
% 10-node quadratic tetrahedron
% 4 Gauss points per element
%
% Input:
%   coords = 10x3 matrix of nodal coordinates
%            rows 1-4  = corner nodes
%            rows 5-10 = midside nodes
%            node order:
%            1,2,3,4   = corners
%            5  = midside 1-2
%            6  = midside 2-3
%            7  = midside 1-3
%            8  = midside 1-4
%            9  = midside 2-4
%            10 = midside 3-4
%
% Output:
%   B_all  = cell(4,1) — one 6x30 B matrix per Gauss point
%   V_all  = 4x1 weights (Ve * wg) for each Gauss point
%   GP_all = 4x3 Gauss point physical coordinates

%% -------------------------------------------------------
% 4-point Gauss quadrature for tetrahedron
% From standard reference tables
%% -------------------------------------------------------

a = 0.1381966011250105;
b = 0.5854101966249685;

% Gauss points in (ξ, η, ζ) reference coords
GP_ref = [a a a;
          b a a;
          a b a;
          a a b];

% Weights (sum = 1/6 for unit tetrahedron)
w = [1/24; 1/24; 1/24; 1/24];

nGP    = 4;
B_all  = cell(nGP, 1);
V_all  = zeros(nGP, 1);
GP_all = zeros(nGP, 3);

%% -------------------------------------------------------
% Loop over Gauss points
%% -------------------------------------------------------

for gp = 1:nGP

    xi  = GP_ref(gp, 1);
    eta = GP_ref(gp, 2);
    zet = GP_ref(gp, 3);

    %% ---------------------------------------------------
    % Tet10 shape functions (quadratic)
    % L1 = 1-xi-eta-zet  (barycentric coord of node 1)
    % L2 = xi
    % L3 = eta
    % L4 = zet
    %% ---------------------------------------------------

    L1 = 1 - xi - eta - zet;
    L2 = xi;
    L3 = eta;
    L4 = zet;

    % Shape functions N1..N10
    % Corner nodes
    % N1  = L1*(2*L1-1)
    % N2  = L2*(2*L2-1)
    % N3  = L3*(2*L3-1)
    % N4  = L4*(2*L4-1)
    % Midside nodes
    % N5  = 4*L1*L2   (edge 1-2)
    % N6  = 4*L2*L3   (edge 2-3)
    % N7  = 4*L1*L3   (edge 1-3)
    % N8  = 4*L1*L4   (edge 1-4)
    % N9  = 4*L2*L4   (edge 2-4)
    % N10 = 4*L3*L4   (edge 3-4)

    %% ---------------------------------------------------
    % Shape function derivatives w.r.t. xi, eta, zet
    % dNi/dxi, dNi/deta, dNi/dzet
    %% ---------------------------------------------------

    % dL1/dxi=-1, dL1/deta=-1, dL1/dzet=-1
    % dL2/dxi= 1, dL2/deta= 0, dL2/dzet= 0
    % dL3/dxi= 0, dL3/deta= 1, dL3/dzet= 0
    % dL4/dxi= 0, dL4/deta= 0, dL4/dzet= 1

    % dNref is 10x3: each row = [dNi/dxi  dNi/deta  dNi/dzet]

    dNref = zeros(10, 3);

    % Corner node 1: N1 = L1*(2L1-1) = (1-xi-eta-zet)*(1-2xi-2eta-2zet)
    % dN1/dxi  = (-1)*(1-2xi-2eta-2zet) + (1-xi-eta-zet)*(-2) = -(4L1-1)
    dNref(1,1) = -(4*L1 - 1);
    dNref(1,2) = -(4*L1 - 1);
    dNref(1,3) = -(4*L1 - 1);

    % Corner node 2: N2 = L2*(2L2-1) = xi*(2xi-1)
    % dN2/dxi = 4xi-1 = 4L2-1
    dNref(2,1) = 4*L2 - 1;
    dNref(2,2) = 0;
    dNref(2,3) = 0;

    % Corner node 3: N3 = L3*(2L3-1) = eta*(2eta-1)
    dNref(3,1) = 0;
    dNref(3,2) = 4*L3 - 1;
    dNref(3,3) = 0;

    % Corner node 4: N4 = L4*(2L4-1) = zet*(2zet-1)
    dNref(4,1) = 0;
    dNref(4,2) = 0;
    dNref(4,3) = 4*L4 - 1;

    % Midside node 5: N5 = 4*L1*L2 = 4*(1-xi-eta-zet)*xi
    % dN5/dxi  = 4*(L1 - L2)  = 4*(1-2xi-eta-zet)
    % dN5/deta = 4*(-xi)       = -4*L2
    % dN5/dzet = 4*(-xi)       = -4*L2
    dNref(5,1) = 4*(L1 - L2);
    dNref(5,2) = -4*L2;
    dNref(5,3) = -4*L2;

    % Midside node 6: N6 = 4*L2*L3 = 4*xi*eta
    % dN6/dxi  = 4*eta = 4*L3
    % dN6/deta = 4*xi  = 4*L2
    % dN6/dzet = 0
    dNref(6,1) = 4*L3;
    dNref(6,2) = 4*L2;
    dNref(6,3) = 0;

    % Midside node 7: N7 = 4*L1*L3 = 4*(1-xi-eta-zet)*eta
    % dN7/dxi  = -4*eta        = -4*L3
    % dN7/deta = 4*(L1 - L3)
    % dN7/dzet = -4*eta        = -4*L3
    dNref(7,1) = -4*L3;
    dNref(7,2) = 4*(L1 - L3);
    dNref(7,3) = -4*L3;

    % Midside node 8: N8 = 4*L1*L4 = 4*(1-xi-eta-zet)*zet
    % dN8/dxi  = -4*zet        = -4*L4
    % dN8/deta = -4*zet        = -4*L4
    % dN8/dzet = 4*(L1 - L4)
    dNref(8,1) = -4*L4;
    dNref(8,2) = -4*L4;
    dNref(8,3) = 4*(L1 - L4);

    % Midside node 9: N9 = 4*L2*L4 = 4*xi*zet
    % dN9/dxi  = 4*zet = 4*L4
    % dN9/deta = 0
    % dN9/dzet = 4*xi  = 4*L2
    dNref(9,1) = 4*L4;
    dNref(9,2) = 0;
    dNref(9,3) = 4*L2;

    % Midside node 10: N10 = 4*L3*L4 = 4*eta*zet
    % dN10/dxi  = 0
    % dN10/deta = 4*zet = 4*L4
    % dN10/dzet = 4*eta = 4*L3
    dNref(10,1) = 0;
    dNref(10,2) = 4*L4;
    dNref(10,3) = 4*L3;

    %% ---------------------------------------------------
    % Jacobian J (3x3)
    % Maps reference coords to physical coords
    % J = dNref' * coords   (3x10 * 10x3 = 3x3)
    %% ---------------------------------------------------

    J    = dNref' * coords;    % 3x3
    detJ = det(J);

    if abs(detJ) < 1e-12
        error('Degenerate Tet10 element — zero Jacobian at GP %d', gp);
    end

    invJ = inv(J);

    %% ---------------------------------------------------
    % Shape function derivatives in physical space
    % dNxyz = dNref * invJ'    (10x3)
    % dNxyz(i,1) = dNi/dx
    % dNxyz(i,2) = dNi/dy
    % dNxyz(i,3) = dNi/dz
    %% ---------------------------------------------------

    dNxyz = dNref * invJ';    % 10x3

    %% ---------------------------------------------------
    % B matrix (6x30)
    % DOF order: [u1 v1 w1  u2 v2 w2 ... u10 v10 w10]
    %% ---------------------------------------------------

    B = zeros(6, 30);

    for i = 1:10

        col = (i-1)*3 + 1;

        dNx = dNxyz(i,1);
        dNy = dNxyz(i,2);
        dNz = dNxyz(i,3);

        B(1, col  ) = dNx;      % exx = du/dx
        B(2, col+1) = dNy;      % eyy = dv/dy
        B(3, col+2) = dNz;      % ezz = dw/dz
        B(4, col  ) = dNy;      % gxy = du/dy + dv/dx
        B(4, col+1) = dNx;
        B(5, col+1) = dNz;      % gyz = dv/dz + dw/dy
        B(5, col+2) = dNy;
        B(6, col  ) = dNz;      % gzx = dw/dx + du/dz
        B(6, col+2) = dNx;

    end

    %% ---------------------------------------------------
    % Physical Gauss point coordinates
    % x_gp = N * coords
    %% ---------------------------------------------------

    N_vec = [L1*(2*L1-1);
             L2*(2*L2-1);
             L3*(2*L3-1);
             L4*(2*L4-1);
             4*L1*L2;
             4*L2*L3;
             4*L1*L3;
             4*L1*L4;
             4*L2*L4;
             4*L3*L4];

    GP_all(gp,:) = N_vec' * coords;

    %% Store
    B_all{gp}  = B;
    V_all(gp)  = abs(detJ) * w(gp);    % weighted volume contribution

end

end