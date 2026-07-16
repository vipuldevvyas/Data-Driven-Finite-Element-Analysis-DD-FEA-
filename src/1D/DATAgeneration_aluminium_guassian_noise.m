% generate_data.m
clear; clc;

% 1. Physical parameters
E_true = 70.0e9;        % True Young's Modulus of Aluminum (70 GPa)
num_points = 20000;     % Density of the data cloud
max_strain = 3.0e-6;    % Upper limit safely above expected max strain (~1.43e-6)
noise_std_pa = 1500;    % Standard deviation of laboratory sensor noise (Pascals)

% 2. Set seed for reproducible random numbers
rng(42);

% 3. Generate random uniform strain coordinates
strain_data = rand(num_points, 1) * max_strain;

% 4. Compute perfect linear-elastic stress and add Gaussian noise
stress_perfect = E_true * strain_data;
experimental_noise = normrnd(0, noise_std_pa, [num_points, 1]);
stress_data = stress_perfect + experimental_noise;

% 5. Create a MATLAB Table and save as CSV
material_table = table(strain_data, stress_data, 'VariableNames', {'Strain', 'Stress_Pa'});
writetable(material_table, 'aluminum_material_data.csv');

fprintf('Success! Material dataset saved as "aluminum_material_data.csv" with %d rows.\n', num_points);

%6. plots

plot(strain_data,stress_data);