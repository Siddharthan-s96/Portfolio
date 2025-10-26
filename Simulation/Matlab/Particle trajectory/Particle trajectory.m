%% Programme: Simplified Winnowing Process Dynamics Program Part 2 
%
% Author: Siddharthan Somasundaram ( matriculation number : 254304 )
%
% Brief description: This program models the dynamics of a simplified 
% winnowing process, which separates grains from chaff based on differences 
% in their mass and aerodynamic properties. Using a combination of Monte-Carlo simulation 
% techniques and Euler's method, the program simulates particle behavior 
% within a controlled airflow environment. Monte-Carlo methods account for 
% stochastic phenomena and variability in particle characteristics, while 
% Euler's method provides deterministic numerical solutions to governing 
% equations of motion. Together, these approaches evaluate the efficiency 
% and dynamics of the winnowing process with a balance of probabilistic and deterministic insights.
%
% Submission date: 07/01/2025
%
%%

clear variables
close all
clc

fprintf('\n***************************************************************************************\n');
fprintf('****************************** Winnowing Process Dynamics *****************************\n');
fprintf('**************************** PROGRAM START time integration ***************************');
fprintf('\n***************************************************************************************\n\n');


%% system parameters

% Initialize invalid particle counters
invalid_grain_count = 0;
invalid_chaff_count = 0;

% Constants and parameters
g = -9.81; % m/s^2 (gravitational acceleration)
rho_fluid = 1.2; % kg/m^3 (density of air)
fluid_mu = 1.8e-5; % Pa*s (dynamic viscosity of air)

% Grain particle properties (values from Table 3.1)
rho_grain = 750; % Density of grain particle
diameter_grain = 2.5e-3; % Grain particle diameter (m)
volume_grain = ((pi/6) * diameter_grain^3);
mass_grain = rho_grain * volume_grain; % Grain particle mass (kg)

% Ensemble Grain particle properties & generation
sd_grain = 1e-3; % Standard deviation of grain diameter
N_grain_particles = 1000; % Number of grain particles
diameter_grain_particles = zeros(N_grain_particles, 1);
mass_grain_particles = zeros(N_grain_particles, 1);

for i = 1:N_grain_particles
    while true
        candidate_diameter = diameter_grain + sd_grain * randn();
        if candidate_diameter > 1e-3
            diameter_grain_particles(i) = candidate_diameter;
            mass_grain_particles(i) = candidate_diameter^3 * rho_grain;
            break;
        else
            invalid_grain_count = invalid_grain_count + 1;
        end
    end
end

% Chaff particle properties
rho_chaff = 50; % Density of chaff particle
diameter_chaff = 3.25e-3; % Chaff particle diameter (m)
volume_chaff = ((pi/6) * diameter_chaff^3);
mass_chaff = rho_chaff * volume_chaff; % Chaff particle mass (kg)

% Ensemble Chaff particle properties & generation
sd_rho_chaff = 20; % Standard deviation of chaff density
min_diameter_chaff = 2e-3; % Minimum chaff diameter (m)
max_diameter_chaff = 5e-3; % Maximum chaff diameter (m)
N_chaff_particles = 1000; % Number of chaff particles

diameter_chaff_particles = zeros(N_chaff_particles, 1);
rho_chaff_particles = zeros(N_chaff_particles, 1);
mass_chaff_particles = zeros(N_chaff_particles, 1);

for i = 1:N_chaff_particles
    while true
        candidate_diameter = min_diameter_chaff + (max_diameter_chaff - min_diameter_chaff) * rand();
        candidate_density = rho_chaff + sd_rho_chaff * randn();
        if candidate_diameter >= min_diameter_chaff && candidate_diameter <= max_diameter_chaff && candidate_density > rho_fluid && candidate_density < rho_grain 
            diameter_chaff_particles(i) = candidate_diameter;
            rho_chaff_particles(i) = candidate_density;
            mass_chaff_particles(i) = candidate_diameter^3 * candidate_density;
            break;
        else
            invalid_chaff_count = invalid_chaff_count + 1;
        end
    end
end

% Simulation parameters
dt = 0.001; % Time step (s)
t_span = 40; % Total simulation time (s)

% Initial conditions
x0 = 0.5; % Initial x-position (m)
y0 = 0.5; % Initial y-position (m)
vx0 = 0; % Initial x-velocity (m/s)
vy0 = 0; % Initial y-velocity (m/s)

% Jet parameters
h = 0.1; % Jet height (m)
u0 = 0.2; % Jet velocity at the inlet (m/s)

% Bin properties
xc = 0.55;
yc = -0.5;
Bin_ground = -0.6;

%number of samples
N_timesamples = 1000;

t1 = tic; % Start the first timer

% Create a waitbar for the simulation of grain and chaff particles
h_wait = waitbar(0, 'Simulating grain particle Trajectory...');

% Simulate grain particle
[t_grain_euler, y_grain_euler] = euler(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, ...
                                        rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);
[t_grain_mc, y_grain_mc] = mc(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, ...
                                        rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground, N_timesamples);

close(h_wait);

h_wait = waitbar(.5, 'Simulating chaff particle Trajectory...');

% Simulate chaff particle
[t_chaff_euler, y_chaff_euler] = euler(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, ...
                                        rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);
[t_chaff_mc, y_chaff_mc] = mc(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, ...
                                        rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground, N_timesamples);

elapsed1 = toc(t1); % Stop the timer and get elapsed time

fprintf('\n***************************************************************************************\n');
fprintf('********************** Trajectory Simulation of Grain and chaff ***********************\n\n');
disp(['Elapsed time for one Grain & Chaff simulation (Euler and Monte-Carlo): ', num2str(elapsed1), ' seconds']);

close(h_wait); % Close the waitbar after completion


% Calculate reference solution with very small dt

t2 = tic; % Start the second timer

% Create a waitbar for the error calculation

dt_ref = 1e-6;

h_wait = waitbar(0, 'Computing Grain reference Trajectory...');

[t_ref_grain, y_grain_ref] = rk4(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, ...
                                    rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt_ref,Bin_ground);

close(h_wait);

h_wait = waitbar(.5, 'Computing Chaff reference Trajectory...');

[t_ref_chaff, y_chaff_ref] = rk4(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, ...
                                        rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt_ref,Bin_ground);

elapsed2 = toc(t2); % Stop the timer and get elapsed time
disp(['Elapsed time for one Grain & Chaff Reference simulation (RK4):', num2str(elapsed2), ' seconds']);

close(h_wait);

%% Calculate trajectories and final positions for 1000 grain and chaff

t3 = tic; % Start the timer
dt = 0.001; % Time step (s)
% Create a waitbar for the calculation trajectories and final positions for 1000 grain

h_wait = waitbar(0, 'Computing Grain Trajectory for 1000 Particles...');

trajectories = cell(N_grain_particles + N_chaff_particles, 1);
final_positions = zeros(N_grain_particles + N_chaff_particles, 3); % [x, y, type (1=grain, 2=chaff)]
theta = zeros(N_grain_particles+N_chaff_particles,1);
terminal_velocity = zeros(N_grain_particles+N_chaff_particles,1);

for i = 1:N_grain_particles
     % Calculate terminal velocity
    [v_terminal,theta_release] = calculate_terminal_velocity(mass_grain_particles(i), diameter_grain_particles(i), ...
                                                             rho_fluid, fluid_mu, g, x0, y0, h, u0);

    terminal_velocity(i) = v_terminal;
    % Generate initial angle for release
    %theta_release = deg2rad(-95 + 10 * rand);

    % Compute initial velocity components
    vx0 = v_terminal * cos(theta_release);
    vy0 = v_terminal * sin(theta_release);

    % Solve trajectory using RK4
    [~, y_traj] = rk4(@(t, y) equations_of_motion(y, mass_grain_particles(i), ...
                diameter_grain_particles(i), rho_fluid, fluid_mu, g, h, u0), ...
                [x0, y0, vx0, vy0], t_span, dt, Bin_ground);

    % Store trajectory and final positions
    trajectories{i} = y_traj;
    final_positions(i, :) = [y_traj(end, 1), y_traj(end, 2), 1]; % Grain particle
    theta(i) = theta_release;
end

close(h_wait); % Close the waitbar after completion
elapsed3 = toc(t3); % Stop the timer and get elapsed time
disp(['Elapsed time for 1000 grain Particle trajectory simulation (RK4): ', num2str(elapsed3), ' seconds']);

t4 = tic; % Start the timer

% Create a waitbar for the calculation trajectories and final positions for 1000 chaff

h_wait = waitbar(.5, 'Computing chaff Trajectory for 1000 Particles...');

for i = 1:N_chaff_particles
    idx = i + N_grain_particles;
    % Calculate terminal velocity
    [v_terminal,theta_release] = calculate_terminal_velocity(mass_chaff_particles(i), diameter_chaff_particles(i), rho_fluid, fluid_mu, g, x0, y0, h, u0);

    terminal_velocity(idx) = v_terminal;
    % Compute initial velocity components
    vx0 = v_terminal * cos(theta_release);
    vy0 = v_terminal * sin(theta_release);

    % Solve trajectory using RK4
    [~, y_traj] = rk4(@(t, y) equations_of_motion(y, mass_chaff_particles(i), ...
                diameter_chaff_particles(i), rho_fluid, fluid_mu, g, h, u0), ...
                [x0, y0, vx0, vy0], t_span, dt, Bin_ground);

    % Store trajectory and final positions
    trajectories{idx} = y_traj;
    final_positions(idx, :) = [y_traj(end, 1), y_traj(end, 2), 2]; % Chaff particle
    theta(idx) = theta_release;
end

close(h_wait); % Close the waitbar after completion
elapsed4 = toc(t4); % Stop the timer and get elapsed time
disp(['Elapsed time for 1000 chaff Particle trajectory simulation (RK4): ', num2str(elapsed4), ' seconds']);


fprintf('---------------------------------------------------------------------------------------\n');
%% Accuracy Gain: Quantify the Improvement of Monte-Carlo Over Euler

fprintf('\n***************************************************************************************');
fprintf('\n***************** Quantifying Accuracy Gain of Monte-Carlo Over Euler *****************\n\n');

t5 = tic; % Start the timer

% Create a waitbar for the calculation trajectories and final positions for 1000 chaff

h_wait = waitbar(.5, 'Calculating Accuracy Gain Over the Trajectory...');

% Interpolate RK4 reference to Euler and MC time steps (Grain)
y_grain_ref_interp_euler = interp1(t_ref_grain, y_grain_ref, t_grain_euler, 'linear', 'extrap');
y_grain_ref_interp_mc = interp1(t_ref_grain, y_grain_ref, t_grain_mc, 'linear', 'extrap');

% Interpolate RK4 reference to Euler and MC time steps (Chaff)
y_chaff_ref_interp_euler = interp1(t_ref_chaff, y_chaff_ref, t_chaff_euler, 'linear', 'extrap');
y_chaff_ref_interp_mc = interp1(t_ref_chaff, y_chaff_ref, t_chaff_mc, 'linear', 'extrap');

% Compute trajectory errors (Grain)
trajectory_error_euler_grain = sqrt(sum((y_grain_ref_interp_euler(:, 1:2) - y_grain_euler(:, 1:2)).^2, 2));
trajectory_error_mc_grain = sqrt(sum((y_grain_ref_interp_mc(:, 1:2) - y_grain_mc(:, 1:2)).^2, 2));

% Compute trajectory errors (Chaff)
trajectory_error_euler_chaff = sqrt(sum((y_chaff_ref_interp_euler(:, 1:2) - y_chaff_euler(:, 1:2)).^2, 2));
trajectory_error_mc_chaff = sqrt(sum((y_chaff_ref_interp_mc(:, 1:2) - y_chaff_mc(:, 1:2)).^2, 2));

% Calculate the point-wise error difference
error_difference_grain = trajectory_error_euler_grain - trajectory_error_mc_grain;
error_difference_chaff = trajectory_error_euler_chaff - trajectory_error_mc_chaff;

% Calculate the *relative* accuracy gain at each time step
relative_accuracy_gain_grain = (error_difference_grain ./ trajectory_error_euler_grain) * 100;
relative_accuracy_gain_chaff = (error_difference_chaff ./ trajectory_error_euler_chaff) * 100;

%Handle cases where Euler error is zero to avoid division by zero
relative_accuracy_gain_grain(trajectory_error_euler_grain == 0) = 0;
relative_accuracy_gain_chaff(trajectory_error_euler_chaff == 0) = 0;

close(h_wait); % Close the waitbar after completion
elapsed5 = toc(t5); % Stop the timer and get elapsed time
disp(['Elapsed time for calculating Accuracy Gain Over the Trajector: ', num2str(elapsed5), ' seconds']);
fprintf('---------------------------------------------------------------------------------------\n');

%% Optimize xc to get grain in Bin 1 and chaff in Bin 2

fprintf('\n***************************************************************************************');
fprintf('\n****************************** Bin Seperator Optimization *****************************\n');

% Initial proportions with xc = 0.55 (showing wrong distributions)
grain_bin2_count_before = sum(final_positions(1:N_grain_particles, 1) > 0.55);
grain_bin2_prop_before = grain_bin2_count_before / N_grain_particles;

chaff_bin1_count_before = sum(final_positions(N_grain_particles+1:end, 1) <= 0.55);
chaff_bin1_prop_before = chaff_bin1_count_before / N_chaff_particles;

% Correct distributions before optimization
%fprintf('\nCorrect distributions before optimization:');
%fprintf(['\nProportion of grain particles in Bin 1 (correct): ', num2str((1 - grain_bin2_prop_before) * 100), '%%']);
%fprintf(['\nProportion of chaff particles in Bin 2 (correct): ', num2str((1 - chaff_bin1_prop_before) * 100), '%%']);

% wrong distributions before optimization
%fprintf('\n\nInitial proportions of particles in wrong bins:');
%fprintf(['\nProportion of grain particles in Bin 2 (should be in Bin 1): ', num2str(grain_bin2_prop_before * 100), '%%']);
%fprintf(['\nProportion of chaff particles in Bin 1 (should be in Bin 2): ', num2str(chaff_bin1_prop_before * 100), '%%']);

t5 = tic; % Start the timer

% Create a waitbar for the finding optimum Xc
h_wait = waitbar(.1, 'Computing optimum Xc for optimium Particle seperation...');

% Optimize xc to get grain in Bin 1 and chaff in Bin 2
xc_values = linspace(min(final_positions(:, 1)), max(final_positions(:, 1)), 5000); % Range of xc values to evaluate
%xc_values = linspace(0, 2, 50); % Range of xc values to evaluate
min_difference = inf;
optimal_xc = xc_values(1);

for xc_temp = xc_values
    % Calculate proportions in wrong bins
    grain_bin2_prop = sum(final_positions(1:N_grain_particles, 1) > xc_temp) / N_grain_particles;
    chaff_bin1_prop = sum(final_positions(N_grain_particles+1:end, 1) <= xc_temp) / N_chaff_particles;
    
    % Find the absolute difference between these proportions
    difference = abs(grain_bin2_prop - chaff_bin1_prop);
    
    if difference < min_difference
        min_difference = difference;
        optimal_xc = xc_temp;
    end
    
end

% Calculate final proportions with optimal xc (showing wrong distributions for comparison)
grain_bin2_count_after = sum(final_positions(1:N_grain_particles, 1) > optimal_xc);
grain_bin2_prop_after = grain_bin2_count_after / N_grain_particles;

chaff_bin1_count_after = sum(final_positions(N_grain_particles+1:end, 1) <= optimal_xc);
chaff_bin1_prop_after = chaff_bin1_count_after / N_chaff_particles;

% display wrong distributions after optimization
%fprintf('\n\nFinal proportions of particles in wrong bins after optimization:');
%fprintf(['\nProportion of grain particles in Bin 2 (should be in Bin 1): ', num2str(grain_bin2_prop_after * 100), '%%']);
%fprintf(['\nProportion of chaff particles in Bin 1 (should be in Bin 2): ', num2str(chaff_bin1_prop_after * 100), '%%']);

% display correct distributions after optimization
%fprintf('\n\nCorrect distributions after optimization:');
%fprintf(['\nProportion of grain particles in Bin 1 (correct): ', num2str((1 - grain_bin2_prop_after) * 100), '%%']);
%fprintf(['\nProportion of chaff particles in Bin 2 (correct): ', num2str((1 - chaff_bin1_prop_after) * 100), '%%']);
%fprintf(['\nOptimal bin boundary (xc) found: ', num2str(optimal_xc)]); 


% Prepare data for the table
description = [
    "Proportion of grain particles in Bin 2 (should be in Bin 1)";
    "Proportion of chaff particles in Bin 1 (should be in Bin 2)";
    "Proportion of grain particles in Bin 1 (correct)";
    "Proportion of chaff particles in Bin 2 (correct)";
];

before_optimization = [
    grain_bin2_prop_before * 100;
    chaff_bin1_prop_before * 100;
    (1 - grain_bin2_prop_before) * 100;
    (1 - chaff_bin1_prop_before) * 100
];

after_optimization = [
    grain_bin2_prop_after * 100;
    chaff_bin1_prop_after * 100;
    (1 - grain_bin2_prop_after) * 100;
    (1 - chaff_bin1_prop_after) * 100
];

% Create a table
result_table = table(description, before_optimization, after_optimization, ...
    'VariableNames', {'Description', 'BeforeOptimization', 'AfterOptimization'});

% Display the table in the command window
fprintf('\nParticle Distribution Proportions Before and After Optimization:\n');
disp(result_table);

% Display the optimal bin boundary
fprintf('\nOptimal bin boundary (Xc) found: %.2f cm\n', optimal_xc);


fprintf('\n---------------------------------------------------------------------------------------\n');

close(h_wait); % Close the waitbar after completion


%% Plot trajectories

% Get the screen size (screen resolution)
screenSize = get(0, 'ScreenSize');

% Define the size of each figure
figWidth = 0.5 * screenSize(3);  % 50% of screen width for the left and right sections
figHeight = (0.5 * screenSize(4)) - 60;  % 50% of screen height for each figure

% Calculate vertical spacing between figures
verticalSpacing = 50;  % Adjust as needed for space between figures


% Define updated colors
grainEulerColor = [0.2, 0.4, 0.8];  % Dark blue
grainMCColor = [0.2, 0.7, 0.2];     % Dark green
chaffEulerColor = [0.6, 0.2, 0.8];  % Purple
chaffMCColor = [0.9, 0.3, 0.2];     % Dark orange-red
binBoundaryColor = [0.9, 0.7, 0.1]; % Gold
jetInletColor = [0.8, 0.2, 0.2];    % Dark red
jetVelocityColor = [0.9, 0.6, 0.1]; % Orange

fig1 = figure(1);

% Trajectories
plot(y_grain_euler(:,1), y_grain_euler(:,2), '--', 'LineWidth', 2, ...
    'Color', grainEulerColor, 'DisplayName', 'Grain Trajectory (Euler)');
hold on;
plot(y_grain_mc(:,1), y_grain_mc(:,2), '-', 'LineWidth', 2, ...
    'Color', grainMCColor, 'DisplayName', 'Grain Trajectory (Monte Carlo)');
plot(y_chaff_euler(:,1), y_chaff_euler(:,2), '--', 'LineWidth', 2, ...
    'Color', chaffEulerColor, 'DisplayName', 'Chaff Trajectory (Euler)');
plot(y_chaff_mc(:,1), y_chaff_mc(:,2), '-', 'LineWidth', 2, ...
    'Color', chaffMCColor, 'DisplayName', 'Chaff Trajectory (Monte Carlo)');

% Bin boundaries (retained values)
x_bin1 = 0.1;       % x-coordinate of the left boundary of Bin 1
x_bin2 = 0.8;       % x-coordinate of the right boundary of Bin 2
x_binsep = xc;      % x-coordinate of the separation between bins
y_bin = yc;         % y-coordinate of the bin separation
y_ground = yc - 0.1; % y-coordinate of the ground

plot([x_bin1, x_bin1, x_binsep, x_binsep, x_binsep, x_binsep, x_bin2, x_bin2], ...
     [y_bin, y_ground, y_ground, y_bin, y_bin, y_ground, y_ground, y_bin], '--', ...
     'LineWidth', 2, 'Color', binBoundaryColor, 'DisplayName', 'Bin Boundaries');

% Wall and jet inlet
plot([0, 0], [-0.6, 0.6], '-', 'LineWidth', 3, 'Color', 'k', 'DisplayName', 'Wall');
plot([0, 0], [-0.05, 0.05], '-', 'LineWidth', 3, ...
    'Color', jetInletColor, 'DisplayName', 'Jet Inlet');

% Bin annotations
text(0.25, -0.55, 'Bin 1', 'FontSize', 14, 'Color', 'b');
text(0.65, -0.55, 'Bin 2', 'FontSize', 14, 'Color', 'b');

% Initial and final particle positions
plot(x0, y0, 'x', 'MarkerSize', 10, 'Color', [0.8500, 0.3250, 0.0980], ...
    'DisplayName', 'Particle Initial Positions');

% Velocity field
x = linspace(0.0001, 0.75, 15);
y = linspace(-0.3, 0.3, 50);
[X, Y] = meshgrid(x, y);
U = 6.2 * u0 * sqrt(h ./ X) .* exp(-50 * (Y.^2) ./ (X.^2));
V = zeros(size(U));
quiver(X, Y, U, V, 1, 'Color', jetVelocityColor, 'DisplayName', 'Jet Velocity');

% Axis limits
xlim([-0.02, 0.81]);
ylim([-0.602, 0.6]);
xlabel('X (m)', 'FontSize', 14);
ylabel('Y (m)', 'FontSize', 14);
title('Particle Trajectory Grain and Chaff', 'FontSize', 16);

% Legend
legend('Location', 'best');
grid on;

% Figure size
set(gcf, 'Position', [figWidth, 50, figWidth, ((figHeight * 2) - 10)]);
hold off;



fig2 = figure(2);

% Define soft, comfortable colors
grainColor = [0.3, 0.6, 0.8];  % Soft teal for grain
chaffColor = [0.6, 0.4, 0.6];  % Soft lavender for chaff

% Subplot for Grain Relative Accuracy Gain
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(t_grain_euler, relative_accuracy_gain_grain, '-', 'LineWidth', 2, 'Color', grainColor);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Relative Accuracy Gain (%)', 'FontSize', 12);
title('Grain: Relative Accuracy Gain of Monte Carlo over Euler', 'FontSize', 14);
grid on;

% Subplot for Chaff Relative Accuracy Gain
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(t_chaff_euler, relative_accuracy_gain_chaff, '-', 'LineWidth', 2, 'Color', chaffColor);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Relative Accuracy Gain (%)', 'FontSize', 12);
title('Chaff: Relative Accuracy Gain of Monte Carlo over Euler', 'FontSize', 14);
grid on;

% Adjust position and size [left, bottom, width, height]
set(gcf, 'Position', [0, 50, figWidth, ((figHeight * 2) - 10)]); 










fig3 = figure(3);

hold on;

title('Trajectories Distribution of 1000 Grain and 1000 Chaff Particles');
xlabel('x position (m)');
ylabel('y position (m)');

% Set axis limits
xlim([-0.02, max(final_positions(:, 1)) + 0.05]);
ylim([-0.62, 0.6]);
set(gcf, 'Position', [0, 50, screenSize(3), screenSize(4)-165]);

% Plot bin boundaries
x_bin1 = 0.1;                                    % x-coordinate of the left boundary of Bin 1
x_bin2 = max(final_positions(:, 1)) + 0.045;    % x-coordinate of the right boundary of Bin 2
x_binsep = xc;                                   % x-coordinate of the separation
y_bin = yc;                                      % y-coordinate of the bin separation
y_ground = yc - 0.102;                           % y-coordinate of the ground

% Boundary of Bin
h_bin_boundary = plot([x_bin1, x_bin1, x_binsep, x_binsep, x_binsep, x_binsep, x_bin2, x_bin2], ...
                      [y_bin, y_ground, y_ground, y_bin, y_bin, y_ground, y_ground, y_bin], '--', ...
                      'LineWidth', 2, 'Color', binBoundaryColor, 'DisplayName', 'Bin Boundaries');

% Boundary of wall
h_wall = plot([0, 0], [-0.6, 0.6], '-', 'LineWidth', 3, 'Color', 'k', 'DisplayName', 'Wall');

% Boundary of slot
h_jet_inlet = plot([0, 0], [-0.05, 0.05], '-', 'LineWidth', 3, ...
                             'Color', jetInletColor, 'DisplayName', 'Jet Inlet');

text(xc / 2, -0.55, 'Bin: 1 ', 'FontSize', 12, 'Color', 'b');
text((xc + max(final_positions(:, 1))) / 2, -0.55, 'Bin: 2 ', 'FontSize', 12, 'Color', 'b');

% Add dynamic text annotations
particle_count_text = text(0.1, 0.4, 'Grain: 0 | Chaff: 0', 'FontSize', 12);
angle_text = text(.1, .35, 'Release Angle: N/A', 'FontSize', 12,'Color', 'b');

% Plot particle positions
Particle_Release = plot(x0, y0, 'o', 'MarkerSize', 2, 'MarkerFaceColor', [0.8500 0.5250 0.0980], ...
    'MarkerEdgeColor', [0.8500 0.5250 0.0980], 'DisplayName', 'Particle Release Position:');

% Define the grid for x and y
x = linspace(0.0001, max(final_positions(:, 1)), max(final_positions(:, 1)) * 20); % x should start from 1 to avoid division by zero
y = linspace(-0.3, 0.3, 50); % y range
[X, Y] = meshgrid(x, y);

% Calculate velocity components
U = 6.2 * u0 * sqrt(h ./ X) .* exp(-50 * (Y .^ 2) ./ (X .^ 2)); % x-component of velocity
V = zeros(size(U)); % y-component of velocity (given as zero in the formula)
jet = quiver(X, Y, U, V, 1, 'Color', jetVelocityColor, 'DisplayName', 'Jet Velocity');

% Initialize counters for particle types
grain_count = 0;
chaff_count = 0;

% Initialize handles for legend entries
grain_handle = []; % For Grain
chaff_handle = []; % For Chaff

% Set separate pause durations for grain and chaff
time_per_particle_grain = 0.00000001; % Pause time for grain particles (in seconds)
time_per_particle_chaff = 0.00000001; % Pause time for chaff particles (in seconds)

for i = 1:length(trajectories)
    set(0, 'CurrentFigure', fig3);
    y_traj = trajectories{i};

    % Skip NaN trajectories
    if any(isnan(y_traj), 'all')
        continue;
    end

    % Extract release angle dynamically from the first time step
    theta_release = atan2d(y_traj(1, 4), y_traj(1, 3)); % Angle in degrees

    % Plot the trajectory
    if final_positions(i, 3) == 1 % Grain
        h_traj = plot(y_traj(:, 1), y_traj(:, 2), 'go', 'MarkerSize', 1, 'LineStyle', 'none');
        grain_count = grain_count + 1; % Increment grain count
        h1 = plot(y_traj(end, 1), y_traj(end, 2), 'g.', 'MarkerSize', 10); % Green for grain

        % Only set DisplayName once for Grain
        if isempty(grain_handle)
            set(h1, 'DisplayName', 'Grain');
            grain_handle = h1;
        end

        % Update dynamic text annotations
        set(particle_count_text, 'String', ['Grain: ', num2str(grain_count), ' | Chaff: ', num2str(chaff_count)]);
        set(angle_text, 'String', ['Angle: ', num2str(theta_release, '%.2f'), '°']);

        % Pause briefly for grain
        pause(time_per_particle_grain);

    else % Chaff
        h_traj = plot(y_traj(:, 1), y_traj(:, 2), 'mo', 'MarkerSize', 1, 'LineStyle', 'none');

        chaff_count = chaff_count + 1; % Increment chaff count
        h2 = plot(y_traj(end, 1), y_traj(end, 2), 'm.', 'MarkerSize', 10); % Magenta for chaff

        % Only set DisplayName once for Chaff
        if isempty(chaff_handle)
            set(h2, 'DisplayName', 'Chaff');
            chaff_handle = h2;
        end

        % Update dynamic text annotations
        set(particle_count_text, 'String', ['Grain: ', num2str(grain_count), ' | Chaff: ', num2str(chaff_count)]);
        set(angle_text, 'String', ['Angle: ', num2str(theta_release, '%.2f'), '°']);

        % Pause briefly for chaff
        pause(time_per_particle_chaff);
    end

    % Delete trajectory after brief display
    delete(h_traj);

end

delete(angle_text);
delete(jet);
% Plot bin boundary
optimal_Bin = line([optimal_xc, optimal_xc], [y_ground, y_bin], 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Optimal Bin Boundary');

legend([grain_handle, chaff_handle, Particle_Release, h_bin_boundary, h_wall, h_jet_inlet, optimal_Bin], 'Location', 'best');


% Set starting positions for the table
x_start = 0.1; % x-coordinate start position for the first column
x_col2 = 0.495;  % x-coordinate start position for the second column
x_col3 = 0.675;  % x-coordinate start position for the third column
y_start = 0.2; % y-coordinate start position
line_height = 0.08; % Vertical spacing between rows

% Title
text(0.5, y_start, 'Proportions Before and After Optimization of Bin boundary Xc', 'FontSize', 14, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
y_start = y_start - line_height; % Move to the next row

% Table Header
text(x_start+.1, y_start, 'Description', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(x_col2, y_start, 'Before Optimization (%)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(x_col3, y_start, 'After Optimization (%)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
y_start = y_start - line_height; % Move to the next row

% Table Rows
for i = 1:length(description)
    % Description column
    text(x_start, y_start, description(i), 'FontSize', 10, 'Interpreter', 'none');
    % Before Optimization column
    text(x_col2, y_start, sprintf('%.1f', before_optimization(i)), 'FontSize', 10, 'HorizontalAlignment', 'center');
    % After Optimization column
    text(x_col3, y_start, sprintf('%.1f', after_optimization(i)), 'FontSize', 10, 'HorizontalAlignment', 'center');
    y_start = y_start - line_height; % Move to the next row
end

text(.1, .35, sprintf('Optimal Bin Boundary: %.2f cm', optimal_xc), 'FontSize', 12, ...
                 'FontWeight', 'bold','Color', [0 0.8 0.5]); 

hold off;

fprintf('\nAll particles simulated.');

fig4 = figure(4);

subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(t_grain_euler, trajectory_error_euler_grain, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Euler Error');
hold on;
plot(t_grain_mc, trajectory_error_mc_grain, 'g-o', 'LineWidth', 1.5, 'DisplayName', 'Monte Carlo Error');
xlabel('Time (s)');
ylabel('Error');
title('Error Comparison vs RK4 Reference for Grain Particle (only for Reference)');
legend('Location', 'best');
grid on;
hold off;

% Subplot for Chaff Error Comparison
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(t_chaff_euler, trajectory_error_euler_chaff, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Euler Error');
hold on;
plot(t_chaff_mc, trajectory_error_mc_chaff, 'g-o', 'LineWidth', 1.5, 'DisplayName', 'Monte Carlo Error');
xlabel('Time (s)');
ylabel('Error');
title('Error Comparison vs RK4 Reference for Chaff Particle (only for Reference)');
legend('Location', 'best');
grid on;
hold off;

set(gcf, 'Position', [figHeight, figWidth/6, figWidth, figHeight*1.5]);

figure(fig1);
figure(fig2);
figure(fig3);

%% Function to calculate drag force
function F_drag = drag_force(velocity_particle_vector, velocity_fluid_vector, diameter, rho_fluid, fluid_mu)

    Cd = Reynolds_number(velocity_particle_vector, velocity_fluid_vector, diameter, rho_fluid, fluid_mu);

    if norm(velocity_fluid_vector - velocity_particle_vector) == 0
         F_drag = zeros(size(velocity_particle_vector));
    else
         F_drag = 0.5 * pi * (diameter^2 / 1) * rho_fluid * Cd * norm(velocity_fluid_vector - velocity_particle_vector) * (velocity_fluid_vector - velocity_particle_vector);
    end
end

%% Function to define jet velocity
function u = jet_velocity(x, y, h, u0)
    if x >= 5 * h
        u = 6.2 * u0 * sqrt(h / x) * exp(-50 * y^2 / x^2);
    else
        u = 0;
    end
end

%% Function to calculate co efficient of drag

function Cd = Reynolds_number(velocity_particle_vector, velocity_fluid_vector, diameter, rho_fluid, fluid_mu)

    Re = rho_fluid * norm(velocity_fluid_vector - velocity_particle_vector) * diameter / fluid_mu; % Reynolds number

    if Re < 800
        Cd = (24 / Re) * (1 + 0.15 * Re^0.687);
    else
        Cd = 0.44;
    end

end

%% Equations of motion
function dydt = equations_of_motion(y, mass, diameter, rho_fluid, fluid_mu, g, h, u0)
    x_pos = y(1);
    y_pos = y(2);
    vx = y(3);
    vy = y(4);

    % Jet velocity
    u_jet = jet_velocity(x_pos, y_pos, h, u0);

    % Relative velocity
    velocity_particle_vector = [vx, vy];
    velocity_fluid_vector = [u_jet, 0];

    % Forces
    F_drag = drag_force(velocity_particle_vector, velocity_fluid_vector, diameter, rho_fluid, fluid_mu);
    F_gravity = [0, mass * g];
    F_buoyancy = [0, -rho_fluid * ((pi / 6) * diameter^3) * g];

    % Accelerations
    net_force = F_gravity + F_buoyancy + F_drag;
    acceleration = net_force / mass;

    dydt = [vx; vy; acceleration(1); acceleration(2)];
   
end

%% Euler method
function [t, y] = euler(odefun, y0, t_span, dt, Bin_ground)
    t = 0:dt:t_span;
    y = zeros(length(t), length(y0));
    y(1, :) = y0;
    for i = 1:length(t) - 1
        dydt = odefun(t(i), y(i, :)');
        y(i + 1, :) = y(i, :) + dydt' * dt;
        if y(i+1, 2) <= Bin_ground % Stop when the particle hits the ground
            break;
        end
    end
    t = t(1:i);
    y = y(1:i, :);
end

%% RK4 method
function [t, y] = rk4(odefun, y0, t_span, dt, Bin_ground)
    t = 0:dt:t_span;
    y = zeros(length(t), length(y0));
    y(1, :) = y0;
    for i = 1:length(t) - 1
        k1 = odefun(t(i), y(i, :)');
        k2 = odefun(t(i) + dt / 2, y(i, :)' + dt * k1 / 2);
        k3 = odefun(t(i) + dt / 2, y(i, :)' + dt * k2 / 2);
        k4 = odefun(t(i) + dt, y(i, :)' + dt * k3);
        dydt = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        y(i + 1, :) = y(i, :) + dydt' * dt;

        if y(i+1, 2) <= Bin_ground
            break;
        end
    end
    t = t(1:i);
    y = y(1:i, :);
end

%% Monte-Carlo Euler Method with Predicted Positions
function [t, y] = mc(odefun, y0, t_span, dt, Bin_ground, N_timesamples)
    % Initialize time and state vectors
    t = 0:dt:t_span;
    y = zeros(length(t), length(y0));
    y(1, :) = y0; % Set initial conditions

    for i = 1:length(t) - 1
        % Step 1: Generate random times within the interval [t_i, t_{i+1}]
        t_samples = t(i) + dt * rand(1, N_timesamples);

        % Step 2: Initialize Monte-Carlo samples
        y_samples = zeros(N_timesamples, length(y0)); % Particle state for each sample

        for j = 1:N_timesamples
            % Step 3: Predict positions using Eq. (3.65)
            x_pred = y(i, 1) + y(i, 3) * (t_samples(j) - t(i)); % x_p = x_p + u_p * dt
            y_pred = y(i, 2) + y(i, 4) * (t_samples(j) - t(i)); % y_p = y_p + v_p * dt

            % Combine predicted positions with current velocities
            y_predicted = [x_pred, y_pred, y(i, 3), y(i, 4)];

            % Step 4: Compute forces using predicted positions
            dydt_predicted = odefun(t_samples(j), y_predicted);

            % Step 5: Use predicted forces to update state
            y_samples(j, :) = y(i, :) + dydt_predicted' * dt;
        end

        % Step 6: Average the sampled trajectories
        y(i + 1, :) = mean(y_samples, 1);

        % Stop if the particle hits the ground
        if y(i + 1, 2) <= Bin_ground
            break;
        end
    end

    % Trim unused entries
    t = t(1:i);
    y = y(1:i, :);
end

%%
function [v_t,theta] = calculate_terminal_velocity(mass, diameter, rho_f, fluid_mu, g, x, y, h, u0)
    % Function to calculate terminal velocity iteratively

    % Initial guess for terminal velocity
    v_t = 1; % Initial guess for terminal velocity magnitude (m/s)
    tol = 1e-6; % Convergence tolerance
    max_iter = 100; % Maximum number of iterations
    theta = deg2rad(-95 + 10 * rand); % Random angle in radians

    for iter = 1:max_iter
        % Random angle for particle velocity vector during the calculation
        
        velocity_particle_vector = [v_t * cos(theta), v_t * sin(theta)]; % Velocity vector

        % Calculate fluid velocity at (x, y) using jet velocity function
        velocity_fluid_vector = [jet_velocity(x, y, h, u0), 0]; % Only x-component for jet

        % Calculate drag coefficient using Reynolds_number
        Cd = Reynolds_number(velocity_particle_vector, velocity_fluid_vector, diameter, rho_f, fluid_mu);

        % Update terminal velocity based on drag and net forces
        v_t_new = sqrt((2 * abs(mass * g)) / (Cd * (pi / 4) * rho_f * diameter^2));

        % Check for convergence
        if abs(v_t_new - v_t) < tol
            break;
        end

        % Update terminal velocity for next iteration
        v_t = v_t_new;
    end
end