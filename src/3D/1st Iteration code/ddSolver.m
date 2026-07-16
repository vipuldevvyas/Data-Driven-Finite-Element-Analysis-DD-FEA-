function results = ddSolver(mesh, K_bc, f_bc, fixedDOFs, ...
                             C, Cinv, DB, F_threshold, maxIter)

% Data-Driven Solver — 3D Linear Elasticity
% Section 3.1, Kirchdoerfer & Ortiz (2015)
% Updated for Tet10 — 4 Gauss points per element

numElem = mesh.numElements;
ndof    = mesh.ndof;
N       = DB.N;
nGP     = numElem * 4;     % 4 GPs per Tet10 element

%% -------------------------------------------------------
% Precompute B matrices and weights for all GPs
%% -------------------------------------------------------

fprintf('\nPrecomputing B matrices...\n');

Bstore  = cell(nGP, 1);
Wstore  = zeros(nGP, 1);
GPdofs  = cell(nGP, 1);

gpCount = 0;

for e = 1:numElem
    conn   = mesh.elements(e,:);
    coords = mesh.nodes(conn,:);
    edofs  = mesh.elementDOFs{e};

    [B_all, V_all, ~] = tet10_Bmatrix(coords);

    for gp = 1:4
        gpCount           = gpCount + 1;
        Bstore{gpCount}   = B_all{gp};
        Wstore(gpCount)   = V_all(gp);
        GPdofs{gpCount}   = edofs;
    end
end

fprintf('Done. Total GPs = %d\n', nGP);

%% -------------------------------------------------------
% Free DOFs
%% -------------------------------------------------------

allDOFs  = 1:ndof;
freeDOFs = setdiff(allDOFs, fixedDOFs);

%% -------------------------------------------------------
% Random initial assignment
%% -------------------------------------------------------

rng(0);
rand_ids = randi(N, nGP, 1);

eps_star = DB.eps(rand_ids, :);
sig_star = DB.sig(rand_ids, :);

%% -------------------------------------------------------
% Storage
%% -------------------------------------------------------

eps_dd   = zeros(nGP, 6);
sig_dd   = zeros(nGP, 6);
Fe       = zeros(nGP, 1);
u_dd     = zeros(ndof, 1);
eta_dd   = zeros(ndof, 1);
Fglobal  = zeros(maxIter, 1);
nChanged = zeros(maxIter, 1);

%% -------------------------------------------------------
% Pre-factorise stiffness matrix once
%% -------------------------------------------------------

fprintf('\nPre-factorising stiffness matrix...\n');
dK = decomposition(K_bc(freeDOFs, freeDOFs), 'lu', 'CheckCondition', false);
fprintf('Done.\n');

fprintf('\n====================================\n');
fprintf('STARTING DD ITERATIONS\n');
fprintf('====================================\n');

%% -------------------------------------------------------
% DD ITERATION LOOP
%% -------------------------------------------------------

for iter = 1:maxIter

    eps_old = eps_star;
    sig_old = sig_star;

    %% ---------------------------------------------------
    % Solve eq. 23a: K * u = Σ Vg * Bᵀ * C * ε*
    %% ---------------------------------------------------

    rhs_u = zeros(ndof, 1);

    for gp = 1:nGP
        B     = Bstore{gp};
        Vg    = Wstore(gp);
        edofs = GPdofs{gp};
        rhs_u(edofs) = rhs_u(edofs) + Vg * B' * C * eps_star(gp,:)';
    end

    u_dd           = zeros(ndof, 1);
    u_dd(freeDOFs) = dK \ rhs_u(freeDOFs);

    %% ---------------------------------------------------
    % Solve eq. 23b: K * η = f - Σ Vg * Bᵀ * σ*
    %% ---------------------------------------------------

    rhs_eta = f_bc;

    for gp = 1:nGP
        B     = Bstore{gp};
        Vg    = Wstore(gp);
        edofs = GPdofs{gp};
        rhs_eta(edofs) = rhs_eta(edofs) - Vg * B' * sig_star(gp,:)';
    end

    eta_dd           = zeros(ndof, 1);
    eta_dd(freeDOFs) = dK \ rhs_eta(freeDOFs);

    %% ---------------------------------------------------
    % Compute local states eq. 20a, 20b
    % ε_gp = B * u
    % σ_gp = σ* + C * B * η
    %% ---------------------------------------------------

    for gp = 1:nGP
        B         = Bstore{gp};
        edofs     = GPdofs{gp};
        u_local   = u_dd(edofs);
        eta_local = eta_dd(edofs);

        eps_dd(gp,:) = (B * u_local)';
        sig_dd(gp,:) = sig_star(gp,:) + (C * B * eta_local)';
    end

    %% ---------------------------------------------------
    % Equilibrium residual
    %% ---------------------------------------------------

    R_int = zeros(ndof, 1);
    for gp = 1:nGP
        B     = Bstore{gp};
        Vg    = Wstore(gp);
        edofs = GPdofs{gp};
        R_int(edofs) = R_int(edofs) + Vg * B' * sig_dd(gp,:)';
    end

    residual      = f_bc - R_int;
    residual_norm = norm(residual(freeDOFs)) / norm(f_bc(freeDOFs));

    %% ---------------------------------------------------
    % Nearest neighbor search — eq. 36 metric
    %% ---------------------------------------------------

    for gp = 1:nGP
        eps_curr = eps_dd(gp,:);
        sig_curr = sig_dd(gp,:);

        dEps = DB.eps - repmat(eps_curr, N, 1);
        dSig = DB.sig - repmat(sig_curr, N, 1);

        dist2 = 0.5*(sum((dEps * C)    .* dEps, 2) + ...
                     sum((dSig * Cinv) .* dSig, 2));

        [bestDist, bestID]  = min(dist2);
        eps_star(gp,:)      = DB.eps(bestID,:);
        sig_star(gp,:)      = DB.sig(bestID,:);
        Fe(gp)              = bestDist;
    end

    %% ---------------------------------------------------
    % Global penalty F = Σ Vg * Fe
    %% ---------------------------------------------------

    Fglobal(iter) = sum(Wstore .* Fe);

    %% ---------------------------------------------------
    % Count changed assignments
    %% ---------------------------------------------------

    changed = 0;
    for gp = 1:nGP
        if any(eps_star(gp,:) ~= eps_old(gp,:)) || ...
           any(sig_star(gp,:) ~= sig_old(gp,:))
            changed = changed + 1;
        end
    end
    nChanged(iter) = changed;

    change_eps = norm(eps_star - eps_old, 'fro') / ...
                 max(norm(eps_star, 'fro'), 1e-20);
    change_sig = norm(sig_star - sig_old, 'fro') / ...
                 max(norm(sig_star, 'fro'), 1e-20);

    fprintf('Iter %3d | F = %.6e | Changed = %4d | eps_rel = %.4e | sig_rel = %.4e | Res = %.4e\n', ...
            iter, Fglobal(iter), nChanged(iter), change_eps, change_sig, residual_norm);

    %% ---------------------------------------------------
    % Convergence checks
    %% ---------------------------------------------------

    if nChanged(iter) == 0
        fprintf('\n====================================\n');
        fprintf('DD-FEM CONVERGED\n');
        fprintf('Iterations = %d\n', iter);
        fprintf('====================================\n');
        break
    end

    if Fglobal(iter) < F_threshold
        fprintf('\nConverged: F below threshold at iter %d\n', iter);
        break
    end

end

%% -------------------------------------------------------
% Store results
%% -------------------------------------------------------

results.u             = u_dd;
results.eta           = eta_dd;
results.eps           = eps_dd;
results.sig           = sig_dd;
results.epsStar       = eps_star;
results.sigStar       = sig_star;
results.Fe            = Fe;
results.Fglobal       = Fglobal(1:iter);
results.nChanged      = nChanged(1:iter);
results.iter          = iter;
results.converged     = (nChanged(iter) == 0);
results.residual_norm = residual_norm;
results.Wstore        = Wstore;
results.nGP           = nGP;

end