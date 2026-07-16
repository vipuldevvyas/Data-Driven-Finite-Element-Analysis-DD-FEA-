function results = ddSolver_KD(mesh, K_bc, f_bc, fixedDOFs, ...
                             C, Cinv, DB, F_threshold, maxIter)

% Data-Driven Solver — 3D Linear Elasticity
% Section 3.1, Kirchdoerfer & Ortiz (2015)
% Updated for Tet10 — 4 Gauss points per element
% MODIFIED: KD-tree nearest-neighbor search (exact, whitened metric)

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
% NEW: Build KD-tree on database using whitened phase-space metric
%
% Theory unchanged — this reproduces eq. (36) of the paper EXACTLY:
%   dist^2 = 0.5*( (eps-eps')*C*(eps-eps')' + (sig-sig')*Cinv*(sig-sig')' )
%
% Whitening: if C = Lc*Lc' and Cinv = Ls*Ls' (Cholesky factors), then
%   ||a*Lc - b*Lc||^2 = (a-b)*C*(a-b)'   exactly
%   ||a*Ls - b*Ls||^2 = (a-b)*Cinv*(a-b)' exactly
% So Euclidean nearest-neighbor in the whitened [eps*Lc, sig*Ls] space
% is IDENTICAL to nearest-neighbor under the original weighted norm.
%% -------------------------------------------------------

fprintf('\nBuilding KD-tree on database (whitened metric)...\n');

Lc = chol(C,    'lower');   % C    = Lc*Lc'
Ls = chol(Cinv, 'lower');   % Cinv = Ls*Ls'

P_db = [DB.eps * Lc, DB.sig * Ls];   % N x 12, whitened database points

Mdl = KDTreeSearcher(P_db, 'Distance', 'euclidean');

fprintf('Done. KD-tree built on %d points.\n', N);

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
    % NEW: Nearest neighbor search via KD-tree — eq. 36 metric
    % Replaces the old per-GP brute-force loop with ONE vectorized
    % batched query for all nGP points at once. Exact same result
    % as brute-force min(dist2), just found via spatial indexing.
    %% ---------------------------------------------------

    P_query = [eps_dd * Lc, sig_dd * Ls];        % nGP x 12, whitened

    [bestID, bestDistWhitened] = knnsearch(Mdl, P_query);   % nGP x 1 each

    eps_star = DB.eps(bestID,:);
    sig_star = DB.sig(bestID,:);

    % Recover the true (unwhitened) penalty value Fe = 0.5 * dist^2
    % knnsearch returns Euclidean distance in whitened space, which by
    % construction equals sqrt(dEps*C*dEps' + dSig*Cinv*dSig'), so:
    Fe = 0.5 * (bestDistWhitened.^2);

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
results.nGP            = nGP;

end