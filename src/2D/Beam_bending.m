clc;
clear;
close all;

fprintf('\n');
fprintf('============================================\n');
fprintf('      2D CANTILEVER FEM BACKBONE\n');
fprintf('============================================\n');
%% =========================================================
% GEOMETRY
%% =========================================================

L = 1.0;        % Beam length (m)
H = 0.1;        % Beam height (m)

nx = 20;
ny = 4;

dx = L/nx;
dy = H/ny;

thickness = 1.0;
%% =========================================================
% MATERIAL
%% =========================================================

E  = 69e9;
nu = 0.33;

C = (E/(1-nu^2))*...
[
    1      nu      0
    nu     1       0
    0      0   (1-nu)/2
];
Cinv = C \ eye(3);
fprintf('\n');
fprintf('Plane Stress Matrix C:\n');
disp(C);
%% =========================================================
% NODE COORDINATES
%% =========================================================

nodes = [];

for j = 0:ny

    for i = 0:nx

        nodes = [nodes;
                 i*dx  j*dy];

    end

end

nNodes = size(nodes,1);

fprintf('\n');
fprintf('Number of Nodes = %d\n',nNodes);
%% =========================================================
% CONNECTIVITY
%% =========================================================

elements = [];

for j = 1:ny

    for i = 1:nx

        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n1 + (nx+1);

        elements = [elements;
                    n1 n2 n3 n4];

    end

end

nElem = size(elements,1);

fprintf('Number of Elements = %d\n',nElem);
%% =========================================================
% GLOBAL DOFs
%% =========================================================

ndof = 2*nNodes;

fprintf('Total DOFs = %d\n',ndof);
%% =========================================================
% GAUSS POINTS
%% =========================================================

g = 1/sqrt(3);

gauss = [
   -g -g
    g -g
    g  g
   -g  g
];
%% =========================================================
% GLOBAL STIFFNESS MATRIX
%% =========================================================

K = sparse(ndof,ndof);
%% =========================================================
% ASSEMBLE GLOBAL K
%% =========================================================

for e = 1:nElem

    elem_nodes = elements(e,:);

    coords = nodes(elem_nodes,:);

    xe = coords(:,1);
    ye = coords(:,2);

    Ke = zeros(8,8);

    for gp = 1:4

        xi  = gauss(gp,1);
        eta = gauss(gp,2);

        dN_dxi = 0.25*[
           -(1-eta)
            (1-eta)
            (1+eta)
           -(1+eta)
        ];

        dN_deta = 0.25*[
           -(1-xi)
           -(1+xi)
            (1+xi)
            (1-xi)
        ];

        J = [
            dN_dxi'*xe    dN_deta'*xe
            dN_dxi'*ye    dN_deta'*ye
        ];

        detJ = det(J);

        invJ = inv(J);

        dN = invJ * [dN_dxi'; dN_deta'];

        dN_dx = dN(1,:);
        dN_dy = dN(2,:);

        B = zeros(3,8);

        for i = 1:4

            B(1,2*i-1) = dN_dx(i);
            B(2,2*i)   = dN_dy(i);

            B(3,2*i-1) = dN_dy(i);
            B(3,2*i)   = dN_dx(i);

        end

        Ke = Ke + ...
            B' * C * B * detJ * thickness;

    end

    dofs = [
        2*elem_nodes(1)-1
        2*elem_nodes(1)
        2*elem_nodes(2)-1
        2*elem_nodes(2)
        2*elem_nodes(3)-1
        2*elem_nodes(3)
        2*elem_nodes(4)-1
        2*elem_nodes(4)
    ];

    K(dofs,dofs) = ...
        K(dofs,dofs) + Ke;

end
%% =========================================================
% LOAD VECTOR
%
% Total load = 1000 N downward
% Distributed equally over the 5 nodes
% on the right edge
%% =========================================================

F = zeros(ndof,1);

%% Right edge nodes

right_nodes = [];

for j = 0:ny

    node = j*(nx+1) + nx + 1;

    right_nodes = [right_nodes; node];

end

fprintf('\n');
fprintf('Right Edge Nodes:\n');
disp(right_nodes);

%% Load per node

P_total = 1000;      % N

P_node = P_total / length(right_nodes);

%% Apply downward force

for i = 1:length(right_nodes)

    node = right_nodes(i);

    %% Vertical DOF
    vdof = 2*node;

    F(vdof) = -P_node;

end

fprintf('\n');
fprintf('Total Applied Force = %.2f N\n',sum(F));

left_nodes = [];

for j = 0:ny

    left_nodes = [left_nodes;
                  j*(nx+1)+1];
end
fixed_dofs = [];

for i = 1:length(left_nodes)

    n = left_nodes(i);

    fixed_dofs = ...
        [fixed_dofs
         2*n-1
         2*n];

end
free_dofs = ...
    setdiff(1:ndof,fixed_dofs);

u = zeros(ndof,1);

u(free_dofs) = ...
    K(free_dofs,free_dofs) \ ...
    F(free_dofs);
%% =========================================================
% TIP DEFLECTION
%% =========================================================

mid_node = right_nodes(3);

tip_deflection = u(2*mid_node);

fprintf('\n');
fprintf('====================================\n');
fprintf('TIP DEFLECTION = %.6e m\n', ...
        tip_deflection);
fprintf('====================================\n');

figure;

scale = 500;

plot(nodes(:,1),...
     nodes(:,2),...
     'k.');

hold on;

deformed = nodes;

deformed(:,1) = ...
    deformed(:,1) + ...
    scale*u(1:2:end);

deformed(:,2) = ...
    deformed(:,2) + ...
    scale*u(2:2:end);

plot(deformed(:,1),...
     deformed(:,2),...
     'ro');

axis equal
grid on

title('Deformed Shape')
%% =========================================================
% DD MATERIAL DATABASE
%% =========================================================

nData = 500000;

strain_scale = 1e-4;

eps_data = strain_scale*(2*rand(nData,3)-1);

sig_data = zeros(nData,3);

for i = 1:nData

    sig_data(i,:) = ...
        (C*eps_data(i,:).')';

end

%% Optional noise

noise_level = 0.0;

sig_data = sig_data .* ...
    (1 + noise_level*randn(size(sig_data)));

fprintf('\n');
fprintf('Database Size = %d\n',nData);
%% =========================================================
% TOTAL GAUSS POINTS
%% =========================================================

nGP = nElem*4;

fprintf('Total GP = %d\n',nGP);
%% =========================================================
% DD STORAGE
%% =========================================================

eps_star = zeros(nGP,3);
sig_star = zeros(nGP,3);

eps_dd = zeros(nGP,3);
sig_dd = zeros(nGP,3);

%% DD objective function

Fe = zeros(nGP,1);


%% =========================================================
% RANDOM INITIAL ASSIGNMENT
%% =========================================================

rand_ids = randi(nData,nGP,1);

eps_star = eps_data(rand_ids,:);

sig_star = sig_data(rand_ids,:);
%% =========================================================
% STORE B MATRICES
%% =========================================================

Bstore = cell(nGP,1);

Wstore = zeros(nGP,1);

GPelem = zeros(nGP,1);

gpCount = 0;

for e = 1:nElem

    elem_nodes = elements(e,:);

    coords = nodes(elem_nodes,:);

    xe = coords(:,1);
    ye = coords(:,2);

    for gp = 1:4

        gpCount = gpCount + 1;

        xi  = gauss(gp,1);
        eta = gauss(gp,2);

        dN_dxi = 0.25*[
           -(1-eta)
            (1-eta)
            (1+eta)
           -(1+eta)
        ];

        dN_deta = 0.25*[
           -(1-xi)
           -(1+xi)
            (1+xi)
            (1-xi)
        ];

        J = [
            dN_dxi'*xe    dN_deta'*xe
            dN_dxi'*ye    dN_deta'*ye
        ];

        detJ = det(J);

        invJ = inv(J);

        dN = invJ * [dN_dxi'; dN_deta'];

        dN_dx = dN(1,:);
        dN_dy = dN(2,:);

        B = zeros(3,8);

        for i = 1:4

            B(1,2*i-1)=dN_dx(i);
            B(2,2*i)=dN_dy(i);

            B(3,2*i-1)=dN_dy(i);
            B(3,2*i)=dN_dx(i);

        end

        Bstore{gpCount}=B;

        Wstore(gpCount)=detJ*thickness;

        GPelem(gpCount)=e;

    end
end
%% =========================================================
% DD REFERENCE STIFFNESS
%% =========================================================

Kdd = K;

%% =========================================================
% DD VARIABLES
%% =========================================================

u_dd = zeros(ndof,1);

eta = zeros(ndof,1);

maxIter = 50;

tol = 1e-10;

%% DD objective storage

Fglobal = zeros(maxIter,1);

nChanged = zeros(maxIter,1);

fprintf('\n');
fprintf('====================================\n');
fprintf('STARTING DD ITERATIONS\n');
fprintf('====================================\n');

for iter = 1:maxIter
    %% =====================================================
    % Store Previous Assignment
    %% =====================================================

    eps_old = eps_star;
    sig_old = sig_star;

    %% =====================================================
    % SECTION 25
    % Solve Equation (11a)
    %
    % K*u = RHS_u
    %% =====================================================

    rhs_u = zeros(ndof,1);

    gpCount = 0;

    for e = 1:nElem

        elem_nodes = elements(e,:);

        dofs = [
            2*elem_nodes(1)-1
            2*elem_nodes(1)
            2*elem_nodes(2)-1
            2*elem_nodes(2)
            2*elem_nodes(3)-1
            2*elem_nodes(3)
            2*elem_nodes(4)-1
            2*elem_nodes(4)
        ];

        for gp = 1:4

            gpCount = gpCount + 1;

            B = Bstore{gpCount};

            rhs_u(dofs) = ...
                rhs_u(dofs) + ...
                Wstore(gpCount) * ...
                B' * C * eps_star(gpCount,:)';

        end

    end

    u_dd = zeros(ndof,1);

    u_dd(free_dofs) = ...
        Kdd(free_dofs,free_dofs) \ ...
        rhs_u(free_dofs);

    %% =====================================================
    % SECTION 26
    % Solve Equation (11b)
    %
    % K*eta = RHS_eta
    %% =====================================================

    rhs_eta = F;

    gpCount = 0;

    for e = 1:nElem

        elem_nodes = elements(e,:);

        dofs = [
            2*elem_nodes(1)-1
            2*elem_nodes(1)
            2*elem_nodes(2)-1
            2*elem_nodes(2)
            2*elem_nodes(3)-1
            2*elem_nodes(3)
            2*elem_nodes(4)-1
            2*elem_nodes(4)
        ];

        for gp = 1:4

            gpCount = gpCount + 1;

            B = Bstore{gpCount};

            rhs_eta(dofs) = ...
                rhs_eta(dofs) - ...
                Wstore(gpCount) * ...
                B' * sig_star(gpCount,:)';

        end

    end

    eta = zeros(ndof,1);

    eta(free_dofs) = ...
        Kdd(free_dofs,free_dofs) \ ...
        rhs_eta(free_dofs);

    %% =====================================================
    % SECTION 27
    % Equation (12)
    % Compute Local States
    %% =====================================================

    gpCount = 0;

    for e = 1:nElem

        elem_nodes = elements(e,:);

        dofs = [
            2*elem_nodes(1)-1
            2*elem_nodes(1)
            2*elem_nodes(2)-1
            2*elem_nodes(2)
            2*elem_nodes(3)-1
            2*elem_nodes(3)
            2*elem_nodes(4)-1
            2*elem_nodes(4)
        ];

        u_local = u_dd(dofs);

        eta_local = eta(dofs);

        for gp = 1:4

            gpCount = gpCount + 1;

            B = Bstore{gpCount};

            eps_dd(gpCount,:) = ...
                (B*u_local)';

            sig_dd(gpCount,:) = ...
                sig_star(gpCount,:) + ...
                (C*(B*eta_local))';

        end

    end

    %% =====================================================
    % SECTION 28
    % Nearest Neighbor Projection
    %% =====================================================

    for gp = 1:nGP

        eps_curr = eps_dd(gp,:)';
        sig_curr = sig_dd(gp,:)';

        bestDist = inf;
        bestID = 1;

        for d = 1:nData

            deps = ...
                eps_curr - ...
                eps_data(d,:)';

            dsig = ...
                sig_curr - ...
                sig_data(d,:)';

            dist = ...
                0.5 * deps' * C * deps + ...
                0.5 * dsig' * Cinv * dsig;

            if dist < bestDist

                bestDist = dist;
                bestID = d;

            end

        end

        eps_star(gp,:) = ...
            eps_data(bestID,:);

        sig_star(gp,:) = ...
            sig_data(bestID,:);

        Fe(gp) = bestDist;
    end

    %% =====================================================
    % SECTION 29
    % Convergence Test
    %% =====================================================
    Fglobal(iter) = ...
    sum(Wstore .* Fe);
    
    %% Relative strain assignment change

change_eps = ...
    norm(eps_star - eps_old,'fro') / ...
    max(norm(eps_star,'fro'),1e-20);

%% Relative stress assignment change

change_sig = ...
    norm(sig_star - sig_old,'fro') / ...
    max(norm(sig_star,'fro'),1e-20);

%% Number of GP assignments changed

changed = 0;

for gp = 1:nGP

    if any(eps_star(gp,:) ~= eps_old(gp,:)) || ...
       any(sig_star(gp,:) ~= sig_old(gp,:))

        changed = changed + 1;

    end

end

nChanged(iter) = changed;

   fprintf( ...
'Iter %3d | F = %.6e | Changed = %4d | eps_rel = %.6e | sig_rel = %.6e\n',...
iter,...
Fglobal(iter),...
nChanged(iter),...
change_eps,...
change_sig);
 
    if nChanged(iter) == 0
        fprintf('\n');
        fprintf('====================================\n');
        fprintf('DD-FEM CONVERGED\n');
        fprintf('Iterations = %d\n',iter);
        fprintf('====================================\n');

        break

    end

end
figure

plot( ...
    1:iter,...
    Fglobal(1:iter),...
    'LineWidth',2);

xlabel('Iteration')
ylabel('F')

title('Global DD Objective Function')

grid on
%% =========================================================
% COMPUTE GAUSS POINT COORDINATES
%% =========================================================

GPcoords = zeros(nGP,2);

gpCount = 0;

for e = 1:nElem

    elem_nodes = elements(e,:);

    coords = nodes(elem_nodes,:);

    xe = coords(:,1);
    ye = coords(:,2);

    for gp = 1:4

        gpCount = gpCount + 1;

        xi  = gauss(gp,1);
        eta = gauss(gp,2);

        N = 0.25 * [
            (1-xi)*(1-eta)
            (1+xi)*(1-eta)
            (1+xi)*(1+eta)
            (1-xi)*(1+eta)
        ];

        xgp = N' * xe;
        ygp = N' * ye;

        GPcoords(gpCount,:) = [xgp ygp];

    end

end
%% =========================================================
% SIGMA XX
%% =========================================================

figure

scatter( ...
    GPcoords(:,1), ...
    GPcoords(:,2), ...
    80, ...
    sig_dd(:,1), ...
    'filled' );

colorbar

xlabel('x (m)')
ylabel('y (m)')

title('\sigma_{xx} (Pa)')

axis equal
grid on
%% =========================================================
% SIGMA YY
%% =========================================================

figure

scatter( ...
    GPcoords(:,1), ...
    GPcoords(:,2), ...
    80, ...
    sig_dd(:,2), ...
    'filled' );

colorbar

xlabel('x (m)')
ylabel('y (m)')

title('\sigma_{yy} (Pa)')

axis equal
grid on
%% =========================================================
% TAU XY
%% =========================================================

figure

scatter( ...
    GPcoords(:,1), ...
    GPcoords(:,2), ...
    80, ...
    sig_dd(:,3), ...
    'filled' );

colorbar

xlabel('x (m)')
ylabel('y (m)')

title('\tau_{xy} (Pa)')

axis equal
grid on
%% =========================================================
% DD-FEM FINAL SUMMARY
%% =========================================================

fprintf('\n');
fprintf('=============================================\n');
fprintf('           DD-FEM FINAL RESULTS\n');
fprintf('=============================================\n');

%% ----------------------------------------------------------
% STRAINS
%% ----------------------------------------------------------

fprintf('\n');
fprintf('STRAINS\n');
fprintf('---------------------------------------------\n');

fprintf('Max eps_xx     = %e\n', max(eps_dd(:,1)));
fprintf('Min eps_xx     = %e\n', min(eps_dd(:,1)));

fprintf('Max eps_yy     = %e\n', max(eps_dd(:,2)));
fprintf('Min eps_yy     = %e\n', min(eps_dd(:,2)));

fprintf('Max gamma_xy   = %e\n', max(eps_dd(:,3)));
fprintf('Min gamma_xy   = %e\n', min(eps_dd(:,3)));

%% ----------------------------------------------------------
% STRESSES
%% ----------------------------------------------------------

fprintf('\n');
fprintf('STRESSES (Pa)\n');
fprintf('---------------------------------------------\n');

fprintf('Max sigma_xx   = %e\n', max(sig_dd(:,1)));
fprintf('Min sigma_xx   = %e\n', min(sig_dd(:,1)));

fprintf('Max sigma_yy   = %e\n', max(sig_dd(:,2)));
fprintf('Min sigma_yy   = %e\n', min(sig_dd(:,2)));

fprintf('Max tau_xy     = %e\n', max(sig_dd(:,3)));
fprintf('Min tau_xy     = %e\n', min(sig_dd(:,3)));

%% ----------------------------------------------------------
% DEFLECTION COMPARISON
%% ----------------------------------------------------------

mid_node = right_nodes(3);

tip_fem = u(2*mid_node);

tip_dd  = u_dd(2*mid_node);

error_percent = ...
    100 * abs(tip_dd - tip_fem) / abs(tip_fem);

fprintf('\n');
fprintf('DEFLECTION COMPARISON\n');
fprintf('---------------------------------------------\n');

fprintf('FEM Tip Deflection      = %e m\n', tip_fem);
fprintf('DD-FEM Tip Deflection   = %e m\n', tip_dd);

fprintf('Absolute Difference     = %e m\n', ...
        abs(tip_dd-tip_fem));

fprintf('Percentage Error        = %.6f %%\n', ...
        error_percent);

%% ----------------------------------------------------------
% MAXIMUM NODAL DEFLECTION
%% ----------------------------------------------------------

v_fem = u(2:2:end);

v_dd  = u_dd(2:2:end);

[maxDef_FEM, idxFEM] = min(v_fem);

[maxDef_DD, idxDD] = min(v_dd);

fprintf('\n');
fprintf('MAXIMUM DEFLECTION\n');
fprintf('---------------------------------------------\n');

fprintf('FEM Max Deflection      = %e m\n', ...
        maxDef_FEM);

fprintf('DD Max Deflection       = %e m\n', ...
        maxDef_DD);

fprintf('Max Deflection Error    = %.6f %%\n', ...
        100*abs(maxDef_DD-maxDef_FEM)/abs(maxDef_FEM));

fprintf('\n');

fprintf('=============================================\n');