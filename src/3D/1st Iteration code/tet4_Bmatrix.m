function [B, V] = tet4_Bmatrix(coords)

% coords = 4x3 matrix of nodal coordinates
% coords(i,:) = [x y z] of node i
%
% Returns:
%   B = 6x12 strain-displacement matrix
%   V = volume of tetrahedron

x = coords(:,1);
y = coords(:,2);
z = coords(:,3);

%% ==========================================================
% JACOBIAN
% Tet4 uses linear shape functions
% J maps reference coords to physical coords
%% ==========================================================

J = [x(2)-x(1)  y(2)-y(1)  z(2)-z(1);
     x(3)-x(1)  y(3)-y(1)  z(3)-z(1);
     x(4)-x(1)  y(4)-y(1)  z(4)-z(1)];

detJ = det(J);

if abs(detJ) < 1e-12
    error('Degenerate tetrahedron — zero volume detected');
end

V = abs(detJ) / 6;

%% ==========================================================
% SHAPE FUNCTION DERIVATIVES
% dN/dx, dN/dy, dN/dz for each of 4 nodes
% These are constant for Tet4 (linear element)
%% ==========================================================

invJ = inv(J);

% Shape function derivatives in reference space
% For Tet4:
% N1 = 1 - xi - eta - zeta
% N2 = xi
% N3 = eta
% N4 = zeta
%
% dN/d(xi,eta,zeta):
dNref = [-1 -1 -1;
          1  0  0;
          0  1  0;
          0  0  1];

% Transform to physical space: dN/d(x,y,z) = dN/d(xi) * inv(J)
dNxyz = dNref * invJ';    % 4x3 matrix
                           % dNxyz(i,1) = dNi/dx
                           % dNxyz(i,2) = dNi/dy
                           % dNxyz(i,3) = dNi/dz

%% ==========================================================
% ASSEMBLE B MATRIX  (6 x 12)
%
% Strain vector: [exx eyy ezz gxy gyz gzx]'
% DOF order per node: [ux uy uz]
% Global DOF order: [u1 v1 w1  u2 v2 w2  u3 v3 w3  u4 v4 w4]
%% ==========================================================

B = zeros(6, 12);

for i = 1:4

    col = (i-1)*3 + 1;    % starting column for node i

    dNx = dNxyz(i,1);
    dNy = dNxyz(i,2);
    dNz = dNxyz(i,3);

    B(1, col  ) = dNx;       % exx = du/dx
    B(2, col+1) = dNy;       % eyy = dv/dy
    B(3, col+2) = dNz;       % ezz = dw/dz
    B(4, col  ) = dNy;       % gxy = du/dy + dv/dx
    B(4, col+1) = dNx;
    B(5, col+1) = dNz;       % gyz = dv/dz + dw/dy
    B(5, col+2) = dNy;
    B(6, col  ) = dNz;       % gzx = dw/dx + du/dz
    B(6, col+2) = dNx;

end

end