clc;
clear;
close all;

%% =========================================================
% TRUE DATA DRIVEN FINITE ELEMENT METHOD (DD-FEM)
%
% Based on:
% Kirchdoerfer-Ortiz DD-FEM formulation
%
% Problem:
% 1D Aluminium Rod
%
% 11 Nodes
% 10 Elements
%
% FEATURES:
% ---------
% 1. Compatibility enforcement
% 2. Equilibrium enforcement
% 3. Local phase-space projection
% 4. Noisy material dataset
% 5. DD-FEM iterative solver
% 6. Stress-strain dataset visualization
% 7. Stress density vs length plot          [NEW]
% 8. Strain density vs length plot          [NEW]
% 9. Deflection vs length (all nodes)       [NEW]
%10. Stress & Strain comparison subplot     [NEW]
%11. Punishment function F vs iterations   [NEW]
%
%% =========================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('      DATA DRIVEN FINITE ELEMENT METHOD\n');
fprintf('============================================\n');

%% =========================================================
% GEOMETRY
%% =========================================================

L  = 1.0;                  % Rod length (m)
A  = 0.01;                 % Cross-sectional area (m^2)

nNodes = 11;
nElem  = 10;

Le = L/nElem;

%% =========================================================
% MATERIAL METRIC PARAMETER
%
% C_e from paper
%% =========================================================

Ce = 70e9;

%% =========================================================
% ELEMENT WEIGHT
%
% w_e = A_e * L_e
%% =========================================================

we = A * Le;

%% =========================================================
% NODAL COORDINATES
%% =========================================================

x = linspace(0,L,nNodes);

%% =========================================================
% CONNECTIVITY MATRIX
%% =========================================================

conn = zeros(nElem,2);

for e = 1:nElem
    conn(e,:) = [e e+1];
end

%% =========================================================
% BUILD GLOBAL STIFFNESS MATRIX
%
% K_ij = sum(w_e*C_e*B_ej*B_ei)
%% =========================================================

K = zeros(nNodes,nNodes);

%% Element B matrix
B = [-1/Le 1/Le];

%% Element stiffness matrix
Ke = we * Ce * (B' * B);

%% Assembly
for e = 1:nElem
    nodes_e = conn(e,:);
    K(nodes_e,nodes_e) = K(nodes_e,nodes_e) + Ke;
end

%% =========================================================
% DISPLAY GLOBAL STIFFNESS MATRIX
%% =========================================================

fprintf('\n');
fprintf('=========== GLOBAL STIFFNESS MATRIX K ===========\n');
disp(K);

%% =========================================================
% EXTERNAL FORCE VECTOR
%% =========================================================

f = zeros(nNodes,1);

%% Applied load at node 11
f(end) = 1000;

%% =========================================================
% BOUNDARY CONDITIONS
%% =========================================================

fixed = 1;
free  = 2:nNodes;

%% =========================================================
% MATERIAL DATASET
%
% Noisy stress-strain cloud
%% =========================================================

nData = 5000;

strain_data = linspace(0,5e-6,nData)';

%% Ideal elastic response
stress_clean = Ce * strain_data;

%% =========================================================
% ADD SMOOTHED NOISE
%% =========================================================

noise_level = 0.02;

noise = noise_level * randn(size(stress_clean));
noise = smoothdata(noise,'gaussian',100);

stress_data = stress_clean .* (1 + noise);

%% Prevent negative stresses
stress_data(stress_data < 0) = 0;

%% =========================================================
% INITIAL LOCAL DATA ASSIGNMENT
%
% Random initialization
%% =========================================================

rand_ids = randi(nData,nElem,1);

eps_star = strain_data(rand_ids);
sig_star = stress_data(rand_ids);

%% =========================================================
% STORAGE VARIABLES
%% =========================================================

u      = zeros(nNodes,1);
eta    = zeros(nNodes,1);
eps_dd = zeros(nElem,1);
sig_dd = zeros(nElem,1);

%% =========================================================
% DD-FEM ITERATION PARAMETERS
%% =========================================================

maxIter = 100;
tol     = 1e-12;

%% =========================================================
% STORAGE FOR PUNISHMENT FUNCTION F vs ITERATIONS  [NEW]
%
% The DD cost (punishment) functional per iteration:
%
%   F(k) = sum_e  w_e * [ (1/2)*C_e*(eps_dd_e - eps_star_e)^2
%                       + (1/2)*(1/C_e)*(sig_dd_e - sig_star_e)^2 ]
%
% This is the total phase-space distance between the
% admissible states and the material dataset states.
%% =========================================================

F_history = zeros(maxIter,1);   % punishment function per iteration
err_history = zeros(maxIter,1); % convergence error per iteration

%% =========================================================
% DD-FEM ITERATIVE SOLVER
%% =========================================================

fprintf('\n');
fprintf('=========== STARTING DD-FEM ITERATION ===========\n');

convergedIter = maxIter;  % default

for k = 1:maxIter

    %% =====================================================
    % STEP 1:
    % Solve displacement problem
    %
    % K*u = rhs_u
    %
    % Equation (11a)
    %% =====================================================

    rhs_u = zeros(nNodes,1);

    for e = 1:nElem
        nodes_e = conn(e,:);
        fe = we * Ce * eps_star(e) * B';
        rhs_u(nodes_e) = rhs_u(nodes_e) + fe;
    end

    %% Apply BC
    Kff     = K(free,free);
    rhs_ff  = rhs_u(free);

    %% Solve for displacement
    u(free)  = Kff \ rhs_ff;
    u(fixed) = 0;

    %% =====================================================
    % STEP 2:
    % Solve Lagrange multiplier system
    %
    % Equation (11b)
    %% =====================================================

    rhs_eta = f;

    for e = 1:nElem
        nodes_e = conn(e,:);
        fe = we * B' * sig_star(e);
        rhs_eta(nodes_e) = rhs_eta(nodes_e) - fe;
    end

    %% Solve for eta
    rhs_eta_f  = rhs_eta(free);
    eta(free)  = Kff \ rhs_eta_f;
    eta(fixed) = 0;

    %% =====================================================
    % STEP 3:
    % Compute admissible local states
    %
    % Equation (12)
    %% =====================================================

    for e = 1:nElem
        nodes_e = conn(e,:);
        ue      = u(nodes_e);
        etae    = eta(nodes_e);
        eps_dd(e) = B * ue;
        sig_dd(e) = sig_star(e) + Ce * (B * etae);
    end

    %% =====================================================
    % STEP 4:
    % Local data projection
    %
    % Find nearest material point
    %% =====================================================

    eps_old = eps_star;
    sig_old = sig_star;

    for e = 1:nElem
        distance = 0.5 * Ce * (eps_dd(e) - strain_data).^2 + ...
                   0.5 * ((sig_dd(e) - stress_data).^2) / Ce;
        [~,idx]      = min(distance);
        eps_star(e)  = strain_data(idx);
        sig_star(e)  = stress_data(idx);
    end

    %% =====================================================
    % STEP 5:
    % Compute Punishment Function F  [NEW]
    %
    % F = sum_e w_e * [ (C_e/2)*(eps_dd_e - eps_star_e)^2
    %                 + (1/(2*C_e))*(sig_dd_e - sig_star_e)^2 ]
    %% =====================================================

    F_val = 0;
    for e = 1:nElem
        F_val = F_val + we * ( ...
            0.5 * Ce  * (eps_dd(e) - eps_star(e))^2 + ...
            0.5 / Ce  * (sig_dd(e) - sig_star(e))^2 );
    end
    F_history(k) = F_val;

    %% =====================================================
    % STEP 6:
    % Convergence check
    %% =====================================================

    err1 = norm(eps_star - eps_old);
    err2 = norm(sig_star - sig_old);
    err  = err1 + err2;

    err_history(k) = err;

    fprintf('Iteration = %3d   Error = %.6e   F = %.6e\n',k,err,F_val);

    if err < tol
        fprintf('\n');
        fprintf('====================================\n');
        fprintf('DD-FEM CONVERGED!\n');
        fprintf('Iterations = %d\n',k);
        fprintf('====================================\n');
        convergedIter = k;
        break;
    end

end

%% Trim history to actual iterations used
F_history   = F_history(1:convergedIter);
err_history = err_history(1:convergedIter);

%% =========================================================
% ELEMENT CENTERS  (must come before analytical fields)
%% =========================================================

xc = zeros(nElem,1);
for e = 1:nElem
    xc(e) = (x(e) + x(e+1))/2;
end

%% =========================================================
% ANALYTICAL SOLUTION
%% =========================================================

P            = 1000;
u_exact      = (P * L) / (A * Ce);
u_analytical = P * x / (A * Ce);

strain_exact = (P/(A*Ce)) * ones(nElem,1);
stress_exact = (P/A)      * ones(nElem,1);

%% =========================================================
% ERRORS
%% =========================================================

stress_error = abs(sig_dd - stress_exact);
strain_error = abs(eps_dd - strain_exact);
disp_error   = abs(u(end) - u_exact);

%% =========================================================
% FINAL RESULTS
%% =========================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('Analytical End Displacement = %.6e m\n', u_exact);
fprintf('DD-FEM End Displacement     = %.6e m\n', u(end));
fprintf('============================================\n');

%% =========================================================
% NODAL / ELEMENT TABLE OUTPUTS
%% =========================================================

nodal_results  = table((1:nNodes)',u,       'VariableNames',{'Node','Displacement_m'});
strain_results = table((1:nElem)', eps_dd,  'VariableNames',{'Element','Strain'});
stress_results = table((1:nElem)', sig_dd,  'VariableNames',{'Element','Stress_Pa'});

fprintf('\n=========== TABLE OUTPUTS ===========\n');
disp(nodal_results);
disp(strain_results);
disp(stress_results);

%% =========================================================
% VALIDATION SUMMARY
%% =========================================================

fprintf('\n========================================\n');
fprintf('VALIDATION SUMMARY\n');
fprintf('========================================\n');
fprintf('Analytical Tip Disp : %.6e m\n',  u_exact);
fprintf('DD-FEM Tip Disp     : %.6e m\n',  u(end));
fprintf('Absolute Error      : %.6e m\n',  disp_error);
fprintf('\n');
fprintf('Mean Stress Error   : %.6e Pa\n', mean(stress_error));
fprintf('Max  Stress Error   : %.6e Pa\n', max(stress_error));
fprintf('\n');
fprintf('Mean Strain Error   : %.6e\n',    mean(strain_error));
fprintf('Max  Strain Error   : %.6e\n',    max(strain_error));
fprintf('========================================\n');

%% =========================================================
%% =========================================================
%%
%%   P L O T   S E C T I O N
%%
%% =========================================================
%% =========================================================

%% Shared plot style
lw = 2.5;
ms = 8;
fs = 12;

%% =========================================================
% PLOT 1  [NEW]
% STRESS DENSITY vs LENGTH
%
% Bar-plot showing stress at each element centroid,
% compared against the analytical uniform stress.
%% =========================================================

figure('Name','Stress Density vs Length','NumberTitle','off');

bar(xc, sig_dd/1e6, 0.6, 'FaceColor',[0.2157 0.4941 0.7216],...
    'EdgeColor','none','FaceAlpha',0.8);

hold on;

plot([0 L],[stress_exact(1) stress_exact(1)]/1e6,...
     'r--','LineWidth',lw,'DisplayName','Analytical');

scatter(xc, sig_dd/1e6, 60, 'ro','filled','DisplayName','DD-FEM');

xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('\sigma (MPa)','FontSize',fs);
title('Stress Density vs Length','FontSize',fs+2,'FontWeight','bold');
legend('DD-FEM (bar)','Analytical','DD-FEM (points)','Location','best');
grid on;
xlim([-0.05 1.05]);

%% =========================================================
% PLOT 2  [NEW]
% STRAIN DENSITY vs LENGTH
%
% Bar-plot showing strain at each element centroid,
% compared against the analytical uniform strain.
%% =========================================================

figure('Name','Strain Density vs Length','NumberTitle','off');

bar(xc, eps_dd, 0.6, 'FaceColor',[0.3020 0.6863 0.2902],...
    'EdgeColor','none','FaceAlpha',0.8);

hold on;

plot([0 L],[strain_exact(1) strain_exact(1)],...
     'r--','LineWidth',lw,'DisplayName','Analytical');

scatter(xc, eps_dd, 60, 'ro','filled','DisplayName','DD-FEM');

xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('\epsilon (dimensionless)','FontSize',fs);
title('Strain Density vs Length','FontSize',fs+2,'FontWeight','bold');
legend('DD-FEM (bar)','Analytical','DD-FEM (points)','Location','best');
grid on;
xlim([-0.05 1.05]);

%% =========================================================
% PLOT 3  [NEW]
% NODAL DEFLECTION vs LENGTH
%
% Displacement at every node (DD-FEM) overlaid with
% the exact linear analytical profile.
%% =========================================================

figure('Name','Nodal Deflection vs Length','NumberTitle','off');

%% Filled area under analytical curve
fill([x fliplr(x)],[u_analytical fliplr(zeros(size(u_analytical)))],...
     [0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.5);

hold on;

plot(x, u_analytical,'k-','LineWidth',lw,'DisplayName','Analytical');

plot(x, u,'b-o','LineWidth',lw,'MarkerSize',ms,...
     'MarkerFaceColor','b','DisplayName','DD-FEM Nodes');

for i = 1:nNodes
    plot([x(i) x(i)],[0 u(i)],'b:','LineWidth',1);
end

xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('Displacement u (m)','FontSize',fs);
title('Nodal Deflection vs Length','FontSize',fs+2,'FontWeight','bold');
legend('Area (Analytical)','Analytical','DD-FEM Nodes','Location','best');
grid on;
xlim([-0.02 1.02]);

%% =========================================================
% PLOT 4  [NEW]
% STRESS AND STRAIN COMPARISON (subplot)
%
% Side-by-side comparison: analytical vs DD-FEM for
% both stress and strain along the rod.
%% =========================================================

figure('Name','Stress & Strain Comparison','NumberTitle','off');

%% --- Stress ---
subplot(2,1,1);

plot(xc, stress_exact/1e6,'k-','LineWidth',lw,'DisplayName','Analytical');
hold on;
plot(xc, sig_dd/1e6,'r-o','LineWidth',lw,'MarkerSize',ms,...
     'MarkerFaceColor','r','DisplayName','DD-FEM');
fill([xc; flipud(xc)],[stress_exact/1e6; flipud(sig_dd/1e6)],...
     [1 0.6 0.6],'EdgeColor','none','FaceAlpha',0.4,'DisplayName','Error band');

xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('\sigma (MPa)','FontSize',fs);
title('Stress: Analytical vs DD-FEM','FontSize',fs,'FontWeight','bold');
legend('Location','best');
grid on;
xlim([0 1]);

%% --- Strain ---
subplot(2,1,2);

plot(xc, strain_exact,'k-','LineWidth',lw,'DisplayName','Analytical');
hold on;
plot(xc, eps_dd,'b-s','LineWidth',lw,'MarkerSize',ms,...
     'MarkerFaceColor','b','DisplayName','DD-FEM');
fill([xc; flipud(xc)],[strain_exact; flipud(eps_dd)],...
     [0.6 0.6 1],'EdgeColor','none','FaceAlpha',0.4,'DisplayName','Error band');

xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('\epsilon (dimensionless)','FontSize',fs);
title('Strain: Analytical vs DD-FEM','FontSize',fs,'FontWeight','bold');
legend('Location','best');
grid on;
xlim([0 1]);

sgtitle('Stress & Strain Comparison Along Rod',...
        'FontSize',fs+3,'FontWeight','bold');

%% =========================================================
% PLOT 5  [NEW]
% PUNISHMENT FUNCTION F vs ITERATIONS
%
%   F(k) = sum_e w_e * [ (C_e/2)*(eps_dd_e - eps_star_e)^2
%                       + (1/(2*C_e))*(sig_dd_e - sig_star_e)^2 ]
%
% Plotted on a semi-log scale to show exponential decay.
%% =========================================================

figure('Name','Punishment Function F vs Iterations','NumberTitle','off');

iterVec = (1:convergedIter)';

%% Semilogy for exponential convergence visibility
semilogy(iterVec, F_history, 'b-o',...
         'LineWidth',lw,'MarkerSize',ms,...
         'MarkerFaceColor','b','DisplayName','F (cost)');

hold on;

%% Overlay convergence error on right y-axis
yyaxis right;
semilogy(iterVec, err_history,'r--s',...
         'LineWidth',lw,'MarkerSize',ms,...
         'MarkerFaceColor','r','DisplayName','||err||');
ylabel('Convergence Error','FontSize',fs,'Color','r');
set(gca,'YColor','r');

yyaxis left;
ylabel('Punishment Function F (J)','FontSize',fs,'Color','b');
set(gca,'YColor','b');

xlabel('Iteration','FontSize',fs);
title('Punishment Function F vs Iterations',...
      'FontSize',fs+2,'FontWeight','bold');

legend('F (cost functional)','Convergence Error','Location','best');

grid on;
xlim([1 convergedIter]);

%% Mark the convergence iteration
if convergedIter < maxIter
    xline(convergedIter,'k--','LineWidth',1.5,...
          'Label',sprintf('Converged at k=%d',convergedIter),...
          'LabelVerticalAlignment','bottom');
end

%% =========================================================
% PLOT 6  (Original)
% DD STATES VS MATERIAL DATASET
%% =========================================================

figure('Name','DD States vs Material Dataset','NumberTitle','off');

scatter(strain_data, stress_data/1e6, 8,'filled',...
        'MarkerFaceAlpha',0.4,'DisplayName','Material Dataset');
hold on;
plot(eps_dd, sig_dd/1e6,'ro','MarkerSize',10,'LineWidth',2,...
     'DisplayName','DD States');
plot(eps_dd, sig_dd/1e6,'r--','LineWidth',1.5,'DisplayName','DD Path');

xlabel('Strain','FontSize',fs);
ylabel('Stress (MPa)','FontSize',fs);
title('DD-FEM States vs Material Dataset','FontSize',fs+2,'FontWeight','bold');
legend('Location','best');
grid on;

%% =========================================================
% PLOT 7  (Original)
% DISPLACEMENT COMPARISON
%% =========================================================

figure('Name','Displacement Comparison','NumberTitle','off');

plot(x, u_analytical,'k-','LineWidth',3,'DisplayName','Analytical');
hold on;
plot(x, u,'ro--','LineWidth',2,'MarkerSize',7,'DisplayName','DD-FEM');

xlabel('Rod Length (m)','FontSize',fs);
ylabel('Displacement (m)','FontSize',fs);
title('Analytical vs DD-FEM Displacement','FontSize',fs+2,'FontWeight','bold');
legend('Location','best');
grid on;

%% =========================================================
% PLOT 8  (Original)
% STRESS COMPARISON
%% =========================================================

figure('Name','Stress Comparison','NumberTitle','off');

plot(xc, stress_exact/1e6,'k-','LineWidth',3,'DisplayName','Analytical');
hold on;
plot(xc, sig_dd/1e6,'ro--','LineWidth',2,'MarkerSize',7,'DisplayName','DD-FEM');

xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('Stress (MPa)','FontSize',fs);
title('Analytical vs DD-FEM Stress','FontSize',fs+2,'FontWeight','bold');
legend('Location','best');
grid on;

%% =========================================================
% PLOT 9  (Original)
% STRAIN COMPARISON
%% =========================================================

figure('Name','Strain Comparison','NumberTitle','off');

plot(xc, strain_exact,'k-','LineWidth',3,'DisplayName','Analytical');
hold on;
plot(xc, eps_dd,'ro--','LineWidth',2,'MarkerSize',7,'DisplayName','DD-FEM');

xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('Strain','FontSize',fs);
title('Analytical vs DD-FEM Strain','FontSize',fs+2,'FontWeight','bold');
legend('Location','best');
grid on;

%% =========================================================
% PLOT 10  (Original)
% STRESS ERROR DISTRIBUTION
%% =========================================================

figure('Name','Stress Error Distribution','NumberTitle','off');
bar(xc, stress_error/1e3,...
    'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none');
xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('Error (kPa)','FontSize',fs);
title('Stress Error Distribution','FontSize',fs+2,'FontWeight','bold');
grid on;

%% =========================================================
% PLOT 11  (Original)
% STRAIN ERROR DISTRIBUTION
%% =========================================================

figure('Name','Strain Error Distribution','NumberTitle','off');
bar(xc, strain_error,...
    'FaceColor',[0.4940 0.1840 0.5560],'EdgeColor','none');
xlabel('Position Along Rod (m)','FontSize',fs);
ylabel('Error','FontSize',fs);
title('Strain Error Distribution','FontSize',fs+2,'FontWeight','bold');
grid on;

%% =========================================================
% END
%% =========================================================
fprintf('\nAll plots generated successfully.\n');