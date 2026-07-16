clear;
clc;
close all;

%% ==========================================================
% INPUT MESH FILE
%% ==========================================================

filename = 'shape1';

%% ==========================================================
% READ FILE
%% ==========================================================

fid = fopen(filename,'r');
if fid == -1
    error('Cannot open mesh file: %s', filename);
end

lines = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
lines = lines{1};

%% ==========================================================
% PHYSICAL GROUPS
%% ==========================================================

physicalNames = containers.Map('KeyType','double','ValueType','char');

idx = find(strcmp(lines,'$PhysicalNames'));

if ~isempty(idx)
    nPhys = str2double(lines{idx+1});
    fprintf('\nPhysical Groups Found:\n');
    for i = 1:nPhys
        txt    = strtrim(lines{idx+1+i});
        tokens = regexp(txt,'(\d+)\s+(\d+)\s+"([^"]+)"','tokens');
        tokens = tokens{1};
        physTag  = str2double(tokens{2});
        physName = tokens{3};
        physicalNames(physTag) = physName;
        fprintf('Tag %d --> %s\n', physTag, physName);
    end
end

%% ==========================================================
% READ NODES
%% ==========================================================

idx = find(strcmp(lines,'$Nodes'));
if isempty(idx)
    idx = find(strcmp(lines,'$ParametricNodes'));
end
if isempty(idx)
    error('No $Nodes or $ParametricNodes section found');
end

numNodes = str2double(lines{idx+1});
nodes    = zeros(numNodes, 3);

for i = 1:numNodes
    vals            = sscanf(lines{idx+1+i}, '%f');
    nodeID          = vals(1);
    nodes(nodeID,:) = vals(2:4)';
end

fprintf('Nodes Read = %d\n', numNodes);

%% ==========================================================
% READ ELEMENTS (MODIFIED: capture Load faces, not just nodes)
%% ==========================================================

idxE = find(strcmp(lines,'$Elements'));
if isempty(idxE)
    error('No $Elements section found');
end

numElementsTotal = str2double(lines{idxE+1});

tetraElements = [];
fixedNodes    = [];
loadNodes     = [];
loadFaces6    = [];   % NEW: Nx6 array, rows = [corner1 corner2 corner3 mid1 mid2 mid3]
loadFacesLin  = [];   % NEW: fallback if only 3-node triangles found on Load group

for i = 1:numElementsTotal

    vals     = sscanf(lines{idxE+1+i}, '%d');
    elemType = vals(2);
    numTags  = vals(3);
    tags     = vals(4 : 3+numTags);
    conn     = vals(4+numTags : end);
    physTag  = tags(1);

    %% Tet10 volume elements
    if elemType == 11
        if length(conn) == 10
            tetraElements = [tetraElements; conn'];
        end

    %% Boundary triangles (linear type 2 OR quadratic type 9)
    elseif elemType == 2 || elemType == 9
        if isKey(physicalNames, physTag)
            groupName = physicalNames(physTag);
            if strcmpi(groupName, 'Fixed')
                fixedNodes = [fixedNodes; conn(:)];
            elseif strcmpi(groupName, 'Load')
                loadNodes  = [loadNodes; conn(:)];
                if elemType == 9 && length(conn) == 6
                    % Gmsh triangle6 order: 1,2,3 corners; 4,5,6 mid-edges
                    loadFaces6 = [loadFaces6; conn'];
                elseif elemType == 2 && length(conn) == 3
                    loadFacesLin = [loadFacesLin; conn'];
                end
            end
        end
    end

end

fixedNodes = unique(fixedNodes);
loadNodes  = unique(loadNodes);

fprintf('Tet10 Elements = %d\n', size(tetraElements,1));
fprintf('Load faces (quadratic, 6-node) found = %d\n', size(loadFaces6,1));
fprintf('Load faces (linear, 3-node) found    = %d\n', size(loadFacesLin,1));

if isempty(loadFaces6) && ~isempty(loadFacesLin)
    warning(['Your mesh only has 3-node (linear) triangles on the Load surface. ' ...
             'Mid-side nodes are NOT available, so consistent Tet10 load distribution ' ...
             'cannot be computed. Fix: in Gmsh, before exporting the .msh file, go to ' ...
             'Tools > Options > Mesh > General, enable "Save all elements", and make sure ' ...
             'the export element order is set to match your volume mesh (2nd order). ' ...
             'Alternatively use File > Export > Mesh with "Element order = 2" and ' ...
             '"Save all elements" checked, so 2D boundary elements are written as 6-node ' ...
             'triangles (elm-type 9), not 3-node (elm-type 2). Falling back to equal split.']);
end

mesh.loadFaces6   = loadFaces6;
mesh.loadFacesLin = loadFacesLin;

%% ==========================================================
% FALLBACK BC DETECTION
%% ==========================================================

if isempty(fixedNodes) || isempty(loadNodes)
    warning('Physical Groups not detected. Using zmin/zmax fallback.');
    tol  = 1e-6;
    zmin = min(nodes(:,3));
    zmax = max(nodes(:,3));
    fixedNodes = find(abs(nodes(:,3) - zmin) < tol);
    loadNodes  = find(abs(nodes(:,3) - zmax) < tol);
end

%% ==========================================================
% GLOBAL DOF NUMBERING
%% ==========================================================

ndof    = 3 * numNodes;
nodeDOF = reshape(1:ndof, 3, numNodes)';

fprintf('Total DOFs = %d\n', ndof);

%% ==========================================================
% ELEMENT DOF TABLE
% Each Tet10 element has 10 nodes x 3 DOFs = 30 DOFs
%% ==========================================================

numElem     = size(tetraElements, 1);
elementDOFs = cell(numElem, 1);

for e = 1:numElem
    conn  = tetraElements(e,:);    % 1x10
    edofs = [];
    for k = 1:10
        edofs = [edofs nodeDOF(conn(k),:)];
    end
    elementDOFs{e} = edofs;        % 1x30
end

%% ==========================================================
% STORE MESH STRUCTURE
%% ==========================================================

mesh.nodes         = nodes;
mesh.elements      = tetraElements;
mesh.numNodes      = numNodes;
mesh.numElements   = numElem;
mesh.fixedNodes    = fixedNodes;
mesh.loadNodes     = loadNodes;
mesh.ndof          = ndof;
mesh.nodeDOF       = nodeDOF;
mesh.elementDOFs   = elementDOFs;
mesh.physicalNames = physicalNames;

%% ==========================================================
% SUMMARY
%% ==========================================================

fprintf('\n====================================\n');
fprintf('MESH SUMMARY\n');
fprintf('====================================\n');
fprintf('Nodes          : %d\n', mesh.numNodes);
fprintf('Tet10 Elements : %d\n', mesh.numElements);
fprintf('Total DOFs     : %d\n', mesh.ndof);
fprintf('Fixed Nodes    : %d\n', length(mesh.fixedNodes));
fprintf('Load Nodes     : %d\n', length(mesh.loadNodes));

%% ==========================================================
% TEST B MATRIX ON FIRST ELEMENT
%% ==========================================================

e      = 1;
conn   = mesh.elements(e,:);
coords = mesh.nodes(conn,:);

[B_all, V_all, GP_all] = tet10_Bmatrix(coords);

fprintf('\nElement 1 — Tet10 check:\n');
fprintf('Total element volume = %.4f mm^3\n', sum(V_all));
fprintf('B matrix size = %d x %d\n', size(B_all{1},1), size(B_all{1},2));
fprintf('Number of Gauss points = %d\n', length(V_all));

%% ==========================================================
% MATERIAL PROPERTIES
%% ==========================================================

fprintf('\n====================================\n');
fprintf('MATERIAL PROPERTIES\n');
fprintf('====================================\n');

E  = input('Enter Young''s Modulus E (MPa): ');
nu = input('Enter Poisson''s Ratio (nu): ');

[C, Cinv] = buildC(E, nu);

lam = E*nu/((1+nu)*(1-2*nu));
mu  = E/(2*(1+nu));

fprintf('\nMaterial Summary\n');
fprintf('------------------------------------\n');
fprintf('Young''s Modulus (E) : %.3f MPa\n', E);
fprintf('Poisson Ratio (nu)  : %.4f\n', nu);
fprintf('Lambda              : %.3f MPa\n', lam);
fprintf('Shear Modulus (mu)  : %.3f MPa\n', mu);
fprintf('\nC matrix (top-left 3x3):\n');
disp(C(1:3,1:3));

%% ==========================================================
% GLOBAL STIFFNESS ASSEMBLY
%% ==========================================================

K = assembleK(mesh, C);

fprintf('K matrix size = %d x %d\n', size(K,1), size(K,2));
fprintf('K is symmetric: %d\n', issymmetric(K));
fprintf('K non-zeros: %d\n', nnz(K));
fprintf('K max diagonal value: %.4e\n', max(full(diag(K))));

%% ==========================================================
% LOAD INPUT (MODIFIED: consistent nodal loads for Tet10 faces)
%% ==========================================================

fprintf('\n====================================\n');
fprintf('LOAD INPUT\n');
fprintf('====================================\n');

Fx = input('Fx (N): ');
Fy = input('Fy (N): ');
Fz = input('Fz (N): ');
F_total = [Fx, Fy, Fz];

f = zeros(ndof, 1);

if ~isempty(mesh.loadFaces6)
    fprintf('Applying CONSISTENT nodal loads (quadratic face integration).\n');

    % --- 1. Compute face areas and total loaded area ---
    numFaces  = size(mesh.loadFaces6, 1);
    faceAreas = zeros(numFaces,1);

    for fidx = 1:numFaces
        corners = mesh.loadFaces6(fidx, 1:3);
        p1 = nodes(corners(1),:);
        p2 = nodes(corners(2),:);
        p3 = nodes(corners(3),:);
        faceAreas(fidx) = 0.5 * norm(cross(p2-p1, p3-p1));
    end
    totalArea = sum(faceAreas);

    fprintf('Total loaded area = %.6f mm^2\n', totalArea);

    % --- 2. Traction vector (uniform, force/area) ---
    traction = F_total / totalArea;

    % --- 3. Distribute via 3-point Gauss quadrature on each face ---
    % Area-coordinate Gauss points (3-pt rule, degree-2 exact)
    gpL = [1/6 1/6 4/6;
           4/6 1/6 1/6;
           1/6 4/6 1/6];
    gpW = [1/3 1/3 1/3];  % weights sum to 1 (will scale by area)

    for fidx = 1:numFaces
        faceConn = mesh.loadFaces6(fidx, :);   % [c1 c2 c3 m1 m2 m3]
        corners  = faceConn(1:3);
        p1 = nodes(corners(1),:);
        p2 = nodes(corners(2),:);
        p3 = nodes(corners(3),:);
        A  = faceAreas(fidx);

        f_local = zeros(6,3);  % local nodal forces [6 nodes x 3 dof]

        for g = 1:3
            L1 = gpL(g,1); L2 = gpL(g,2); L3 = gpL(g,3);

            % Quadratic triangle shape functions
            N = zeros(6,1);
            N(1) = L1*(2*L1-1);
            N(2) = L2*(2*L2-1);
            N(3) = L3*(2*L3-1);
            N(4) = 4*L1*L2;
            N(5) = 4*L2*L3;
            N(6) = 4*L3*L1;

            w = gpW(g) * A;  % weight * area (Jacobian already folded into A for flat tri)

            for k = 1:6
                f_local(k,:) = f_local(k,:) + N(k) * traction * w;
            end
        end

        % Assemble into global force vector
        for k = 1:6
            nodeID = faceConn(k);
            f(nodeDOF(nodeID,1)) = f(nodeDOF(nodeID,1)) + f_local(k,1);
            f(nodeDOF(nodeID,2)) = f(nodeDOF(nodeID,2)) + f_local(k,2);
            f(nodeDOF(nodeID,3)) = f(nodeDOF(nodeID,3)) + f_local(k,3);
        end
    end

    fprintf('\nSanity check — corner vs mid-side node force split:\n');
    cornerIDs = unique(mesh.loadFaces6(:,1:3));
    midIDs    = unique(mesh.loadFaces6(:,4:6));
    fprintf('Sum |Fx| on corner nodes  : %.6e N\n', sum(abs(f(nodeDOF(cornerIDs,1)))));
    fprintf('Sum |Fx| on mid-side nodes: %.6e N\n', sum(abs(f(nodeDOF(midIDs,1)))));
    fprintf('(Expect corner ~0, mid-side ~Fx for flat/uniform load)\n');

else
    warning('No quadratic (6-node) load faces available — falling back to equal split across all load nodes.');
    numLoadNodes = length(mesh.loadNodes);
    for i = 1:numLoadNodes
        nodei = mesh.loadNodes(i);
        f(nodeDOF(nodei,1)) = f(nodeDOF(nodei,1)) + Fx/numLoadNodes;
        f(nodeDOF(nodei,2)) = f(nodeDOF(nodei,2)) + Fy/numLoadNodes;
        f(nodeDOF(nodei,3)) = f(nodeDOF(nodei,3)) + Fz/numLoadNodes;
    end
end

fprintf('\nTotal Applied Force\n');
fprintf('Fx = %.4f N\n', sum(f(1:3:end)));
fprintf('Fy = %.4f N\n', sum(f(2:3:end)));
fprintf('Fz = %.4f N\n', sum(f(3:3:end)));

%% ==========================================================
% BOUNDARY CONDITIONS
%% ==========================================================

K_bc = K;
f_bc = f;

fixedDOFs = [];
for i = 1:length(mesh.fixedNodes)
    nodei     = mesh.fixedNodes(i);
    fixedDOFs = [fixedDOFs mesh.nodeDOF(nodei,:)];
end
fixedDOFs = unique(fixedDOFs);

for i = 1:length(fixedDOFs)
    dof = fixedDOFs(i);
    K_bc(dof,:)    = 0;
    K_bc(:,dof)    = 0;
    K_bc(dof,dof)  = 1;
    f_bc(dof)      = 0;
end

fprintf('\nFixed DOFs: %d\n', length(fixedDOFs));
fprintf('K_bc symmetric: %d\n', issymmetric(K_bc));

%% ==========================================================
% CLASSICAL FEM SOLVE
%% ==========================================================

u = K_bc \ f_bc;

fprintf('\nDisplacement solution:\n');
fprintf('Max displacement : %.6e mm\n', max(abs(u)));
fprintf('Max X disp       : %.6e mm\n', max(abs(u(1:3:end))));
fprintf('Max Y disp       : %.6e mm\n', max(abs(u(2:3:end))));
fprintf('Max Z disp       : %.6e mm\n', max(abs(u(3:3:end))));

U = reshape(u, 3, mesh.numNodes)';

fprintf('\nSample displacements (first 5 nodes):\n');
fprintf('Node   Ux            Uy            Uz\n');
for i = 1:5
    fprintf('%4d   %12.4e  %12.4e  %12.4e\n', i, U(i,1), U(i,2), U(i,3));
end

mesh.u = u;
mesh.U = U;

%% ==========================================================
% STRESS AND STRAIN RECOVERY — Tet10
% Loop over elements, loop over 4 GPs per element
%% ==========================================================

nGP_total = numElem * 4;

strain   = zeros(nGP_total, 6);
stress   = zeros(nGP_total, 6);
GPcoords = zeros(nGP_total, 3);
GPvols   = zeros(nGP_total, 1);
GPelems  = zeros(nGP_total, 1);

gpCount = 0;

for e = 1:numElem

    conn   = mesh.elements(e,:);
    coords = mesh.nodes(conn,:);
    edofs  = mesh.elementDOFs{e};
    ue     = u(edofs);

    [B_all, V_all, GP_all] = tet10_Bmatrix(coords);

    for gp = 1:4
        gpCount = gpCount + 1;
        B       = B_all{gp};

        eps_gp  = B * ue;
        sig_gp  = C * eps_gp;

        strain(gpCount,:)   = eps_gp';
        stress(gpCount,:)   = sig_gp';
        GPcoords(gpCount,:) = GP_all(gp,:);
        GPvols(gpCount)     = V_all(gp);
        GPelems(gpCount)    = e;
    end

end

mesh.strain   = strain;
mesh.stress   = stress;
mesh.GPcoords = GPcoords;
mesh.GPvols   = GPvols;

fprintf('\n====================================\n');
fprintf('STRESS/STRAIN SUMMARY (Tet10)\n');
fprintf('====================================\n');

fprintf('\nMax strain components:\n');
fprintf('exx: %12.4e\n', max(abs(strain(:,1))));
fprintf('eyy: %12.4e\n', max(abs(strain(:,2))));
fprintf('ezz: %12.4e\n', max(abs(strain(:,3))));
fprintf('gxy: %12.4e\n', max(abs(strain(:,4))));
fprintf('gyz: %12.4e\n', max(abs(strain(:,5))));
fprintf('gzx: %12.4e\n', max(abs(strain(:,6))));

fprintf('\nMax stress components (MPa):\n');
fprintf('sxx: %12.4e\n', max(abs(stress(:,1))));
fprintf('syy: %12.4e\n', max(abs(stress(:,2))));
fprintf('szz: %12.4e\n', max(abs(stress(:,3))));
fprintf('txy: %12.4e\n', max(abs(stress(:,4))));
fprintf('tyz: %12.4e\n', max(abs(stress(:,5))));
fprintf('tzx: %12.4e\n', max(abs(stress(:,6))));

%% Nodal stress averaging
nodalStress = zeros(mesh.numNodes, 6);
nodalCount  = zeros(mesh.numNodes, 1);

for e = 1:numElem
    conn    = mesh.elements(e,:);
    gpStart = (e-1)*4 + 1;
    gpEnd   = e*4;
    avgStress = mean(stress(gpStart:gpEnd,:), 1);
    for k = 1:10
        nodalStress(conn(k),:) = nodalStress(conn(k),:) + avgStress;
        nodalCount(conn(k))    = nodalCount(conn(k)) + 1;
    end
end

for i = 1:mesh.numNodes
    if nodalCount(i) > 0
        nodalStress(i,:) = nodalStress(i,:) / nodalCount(i);
    end
end

mesh.nodalStress = nodalStress;

fprintf('\nNodal stress (smoothed) max values:\n');
fprintf('sxx: %.4f MPa\n', max(nodalStress(:,1)));
fprintf('syy: %.4f MPa\n', max(nodalStress(:,2)));
fprintf('szz: %.4f MPa\n', max(nodalStress(:,3)));

%% Von Mises at Gauss points
s = stress;
vonMises = sqrt(0.5*((s(:,1)-s(:,2)).^2 + ...
                      (s(:,2)-s(:,3)).^2 + ...
                      (s(:,3)-s(:,1)).^2 + ...
                      6*(s(:,4).^2+s(:,5).^2+s(:,6).^2)));

mesh.vonMises = vonMises;
fprintf('\nVon Mises (GP):\n');
fprintf('Max : %12.4e MPa\n', max(vonMises));
fprintf('Min : %12.4e MPa\n', min(vonMises));
fprintf('Mean: %12.4e MPa\n', mean(vonMises));

%% ==========================================================
% DATABASE SETTINGS
%% ==========================================================

fprintf('\n====================================\n');
fprintf('DATABASE SETTINGS\n');
fprintf('====================================\n');

N          = input('Enter Number of Database Points : ');
noiseLevel = 0.0;

fprintf('\nDatabase Size = %d\n', N);
fprintf('Noise Level   = %.2f\n', noiseLevel);

DB = generateDatabase(C, Cinv, strain, N, noiseLevel);

%% ==========================================================
% DD SOLVER
%% ==========================================================

F_threshold = 1e-6;
maxIter     = 200;

results = ddSolver(mesh, K_bc, f_bc, fixedDOFs, ...
                   C, Cinv, DB, F_threshold, maxIter);

%% ==========================================================
% DD RESULTS SUMMARY
%% ==========================================================

fprintf('\n=============================================\n');
fprintf('           DD-FEM FINAL RESULTS\n');
fprintf('=============================================\n');

fprintf('\nSTRAINS\n');
fprintf('---------------------------------------------\n');
fprintf('Max exx  = %e\n', max(results.eps(:,1)));
fprintf('Min exx  = %e\n', min(results.eps(:,1)));
fprintf('Max eyy  = %e\n', max(results.eps(:,2)));
fprintf('Min eyy  = %e\n', min(results.eps(:,2)));
fprintf('Max ezz  = %e\n', max(results.eps(:,3)));
fprintf('Min ezz  = %e\n', min(results.eps(:,3)));
fprintf('Max gxy  = %e\n', max(results.eps(:,4)));
fprintf('Max gyz  = %e\n', max(results.eps(:,5)));
fprintf('Max gzx  = %e\n', max(results.eps(:,6)));

fprintf('\nSTRESSES (MPa)\n');
fprintf('---------------------------------------------\n');
fprintf('Max sxx  = %e\n', max(results.sig(:,1)));
fprintf('Min sxx  = %e\n', min(results.sig(:,1)));
fprintf('Max syy  = %e\n', max(results.sig(:,2)));
fprintf('Min syy  = %e\n', min(results.sig(:,2)));
fprintf('Max szz  = %e\n', max(results.sig(:,3)));
fprintf('Min szz  = %e\n', min(results.sig(:,3)));

fprintf('\nDEFLECTION COMPARISON\n');
fprintf('---------------------------------------------\n');
fprintf('FEM Max disp    = %e mm\n', max(abs(u)));
fprintf('DD  Max disp    = %e mm\n', max(abs(results.u)));
err = 100*abs(max(abs(results.u)) - max(abs(u))) / max(abs(u));
fprintf('Error           = %.4f %%\n', err);

fprintf('\nDD Residual norm = %.4e\n', results.residual_norm);


%% ==========================================================
% PLOTS
%% ==========================================================

%% Convergence plot
figure;
plot(1:results.iter, results.Fglobal, 'b-o', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Iteration');
ylabel('F (Global Penalty)');
title('DD Global Objective Function (eq. 17)');
grid on;

%% Mesh nodes
figure;
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 15, 'k', 'filled');
hold on;
scatter3(nodes(fixedNodes,1), nodes(fixedNodes,2), nodes(fixedNodes,3), 80, 'b', 'filled');
scatter3(nodes(loadNodes,1),  nodes(loadNodes,2),  nodes(loadNodes,3),  80, 'r', 'filled');
grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Mesh Nodes');
legend('All Nodes','Fixed Nodes','Load Nodes');

%% Stress contour at Gauss points — sxx
figure;
scatter3(GPcoords(:,1), GPcoords(:,2), GPcoords(:,3), ...
         30, stress(:,1), 'filled');
colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('\sigma_{xx} at Gauss Points (MPa)');
axis equal; grid on;

%% Stress contour — syy
figure;
scatter3(GPcoords(:,1), GPcoords(:,2), GPcoords(:,3), ...
         30, stress(:,2), 'filled');
colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('\sigma_{yy} at Gauss Points (MPa)');
axis equal; grid on;

%% Stress contour — szz
figure;
scatter3(GPcoords(:,1), GPcoords(:,2), GPcoords(:,3), ...
         30, stress(:,3), 'filled');
colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('\sigma_{zz} at Gauss Points (MPa)');
axis equal; grid on;

%% Von Mises contour
figure;
scatter3(GPcoords(:,1), GPcoords(:,2), GPcoords(:,3), ...
         30, vonMises, 'filled');
colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Von Mises Stress at Gauss Points (MPa)');
axis equal; grid on;