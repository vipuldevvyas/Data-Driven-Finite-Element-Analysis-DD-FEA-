function DB = generateDatabase(C, Cinv, strain_ref, N, noiseLevel)

% Generates synthetic material database
% Mirrors the 2D reference code database generation
% Points sampled near constitutive law + Gaussian noise
%
% For 3D: strain/stress vectors are 6-component

%% Determine strain range from reference FEM solution
epsMax = max(abs(strain_ref(:))) * 1.5;

if epsMax < 1e-12
    epsMax = 1e-6;
end

fprintf('\n--- Database Generation ---\n');
fprintf(['Strain r' ...
    'ange  : +/- %.4e\n'], epsMax);
fprintf('Database size : %d points\n', N);
fprintf('Noise level   : %.4f\n', noiseLevel);

%% Sample strain points uniformly in 6D
rng(42);

eps_data = epsMax * (2*rand(N,6) - 1);    % N x 6

%% Compute exact stress from constitutive law
% sig = C * eps  (on the constitutive manifold)
sig_data = (C * eps_data')';              % N x 6

%% Add Gaussian noise to both strain and stress (additive, matching paper §2.2 / Fig. 7)
% Paper adds noise additively with variance ρ_k, not multiplicatively.
% Noise is scaled to the range of each field so the magnitude is meaningful.
epsMax_field = max(abs(eps_data(:)));
sigMax       = max(abs(sig_data(:)));

eps_data = eps_data + noiseLevel * epsMax_field * randn(size(eps_data));
sig_data = sig_data + noiseLevel * sigMax       * randn(size(sig_data));

%% Store
DB.eps    = eps_data;
DB.sig    = sig_data;
DB.N      = N;
DB.epsMax = epsMax;
DB.sigMax = sigMax;

fprintf('Stress range  : +/- %.4e MPa\n', sigMax);
fprintf('Database ready.\n');

end