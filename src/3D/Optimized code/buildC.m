function [C, Cinv] = buildC(E, nu)

% Builds 6x6 isotropic linear elastic material matrix
% and its inverse for the DD distance metric
%
% Inputs:
%   E   = Young's modulus (Pa or MPa, consistent with your units)
%   nu  = Poisson's ratio
%
% Outputs:
%   C    = 6x6 stiffness matrix  (stress = C * strain)
%   Cinv = 6x6 compliance matrix (strain = Cinv * stress)

lam = E*nu / ((1+nu)*(1-2*nu));     % Lame lambda
mu  = E / (2*(1+nu));                % Shear modulus

C = [lam+2*mu  lam      lam      0   0   0;
     lam        lam+2*mu lam      0   0   0;
     lam        lam      lam+2*mu 0   0   0;
     0          0        0        mu  0   0;
     0          0        0        0   mu  0;
     0          0        0        0   0   mu];

Cinv = inv(C);

end