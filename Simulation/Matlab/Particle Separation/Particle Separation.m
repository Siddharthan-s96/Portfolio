
%% Programme: Simplified Winnowing Process Dynamics Program 
%
% Author: Siddharthan Somasundaram ( matriculation number : 254304 )
%
% Brief description: This program models the dynamics of a simplified 
% winnowing process, designed to separate grains from chaff based on their
% differences in mass and aerodynamic properties, By simulating the
% behavior of particles in a controlled airflow.
%
% Submission date: 25/11/2024
%
%%

clear variables
close all
clc

fprintf('\n*********************************************************\n');
fprintf('*************** Winnowing Process Dynamics **************\n');
fprintf('************* PROGRAM START time integration ************');
fprintf('\n*********************************************************\n\n');

%% system parameters

% Constants and parameters
g = -9.81; % m/s^2 (gravitational acceleration)
rho_fluid = 1.2; % kg/m^3 (density of air)
fluid_mu = 1.8e-5; % Pa*s (dynamic viscosity of air)

% Grain particle properties (values from Table 3.1)
rho_grain = 750; % Density of grain particle
diameter_grain = 2.5e-3; % Grain particle diameter (m)
volume_grain = ((pi/6) * diameter_grain^3);
mass_grain = rho_grain * volume_grain; % Grain particle mass (kg)

% Chaff particle properties
rho_chaff = 50; % Density of chaff particle
diameter_chaff = 3.25e-3; % Chaff particle diameter (m)
volume_chaff = ((pi/6) * diameter_chaff^3);
mass_chaff = rho_chaff * volume_chaff; % Chaff particle mass (kg)

% Simulation parameters
dt = 0.01; % Time step (s)
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

t1 = tic; % Start the first timer

% Create a waitbar for the simulation of grain and chaff particles
h_wait = waitbar(0, 'Simulating grain particle Trajectory...');

% Simulate grain particle
[t_grain_euler, y_grain_euler] = euler(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);
[t_grain_rk4, y_grain_rk4] = rk4(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);

close(h_wait);

h_wait = waitbar(50, 'Simulating chaff particle Trajectory...');

% Simulate chaff particle
[t_chaff_euler, y_chaff_euler] = euler(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);
[t_chaff_rk4, y_chaff_rk4] = rk4(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);

elapsed1 = toc(t1); % Stop the timer and get elapsed time
disp(['Elapsed time for Particle simulation: ', num2str(elapsed1), ' seconds']);

close(h_wait); % Close the waitbar after completion

%Error Calculation
% Calculate reference solution with very small dt

t2 = tic; % Start the second timer

% Create a waitbar for the error calculation

dt_values = [.01,.001,.0001,.00001]; % Different time steps to test
dt_ref = 1e-6;

h_wait = waitbar(0, 'Computing Grain Actual Trajectory...');

[t_ref_grain, y_grain_ref] = rk4(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt_ref,Bin_ground);

close(h_wait);

h_wait = waitbar(.125, 'Computing Chaff Actual Trajectory...');

[t_ref_chaff, y_chaff_ref] = rk4(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt_ref,Bin_ground);

close(h_wait);

% Preallocate error arrays
error_euler_grain = zeros(length(dt_values), 1);
error_rk4_grain = zeros(length(dt_values), 1);
error_euler_chaff = zeros(length(dt_values), 1);
error_rk4_chaff = zeros(length(dt_values), 1);

for i = 1:length(dt_values)
    dt = dt_values(i);
    
    h_wait = waitbar(i * .25, 'Calculating error for different dt values...');

    % Simulate heavy particle using Euler method
    [t_grain_euler_error, y_grain_euler_error] = euler(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);
   
    % Simulate heavy particle using RK4 method
    [t_grain_rk4_error, y_grain_rk4_error] = rk4(@(t, y) equations_of_motion(y, mass_grain, diameter_grain, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt,Bin_ground);
   
    % Simulate light particle using Euler method
    [t_chaff_euler_error, y_chaff_euler_error] = euler(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);
   
    % Simulate light particle using RK4 method
    [t_chaff_rk4_error, y_chaff_rk4_error] = rk4(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);

    % Interpolate reference solutions to match time points
    y_ref_grain_euler = interp1(t_ref_grain, y_grain_ref, t_grain_euler_error);
    y_ref_grain_rk4 = interp1(t_ref_grain, y_grain_ref, t_grain_rk4_error);
    y_ref_chaff_euler = interp1(t_ref_chaff, y_chaff_ref, t_chaff_euler_error);
    y_ref_chaff_rk4 = interp1(t_ref_chaff, y_chaff_ref, t_chaff_rk4_error);

    % Compute relative error using the provided formula
    error_euler_grain(i) = norm(y_grain_euler_error(end, 1:2) - y_ref_grain_euler(end, 1:2));    % Relative error for Euler grain
    error_rk4_grain(i) = norm(y_grain_rk4_error(end, 1:2) - y_ref_grain_rk4(end, 1:2));          % Relative error for RK4 grain

    error_euler_chaff(i) = norm(y_chaff_euler_error(end, 1:2) - y_ref_chaff_euler(end, 1:2));    % Relative error for Euler chaff
    error_rk4_chaff(i) = norm(y_chaff_rk4_error(end, 1:2) - y_ref_chaff_rk4(end, 1:2));          % Relative error for RK4 chaff

    close(h_wait); % Close the waitbar after completion
end

elapsed2 = toc(t2); % Stop the timer and get elapsed time
disp(['Elapsed time for Error Calculation: ', num2str(elapsed2), ' seconds']);

%Minimum Jet Velocity for Particle seperation

t3 = tic; % Start the second timer

optimal_u0 = mim_jet_velocity(mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0, x0, y0, vx0, vy0, xc, yc, t_span, dt, Bin_ground);
[t_chaff_optimal_rk4, y_chaff_optimal_rk4] = rk4(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, optimal_u0), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);

elapsed3 = toc(t3); % Stop the timer and get elapsed time
disp(['Elapsed time to find Optimal(minimum) Velocity: ', num2str(elapsed3), ' seconds']);
fprintf('Optimal initial jet velocity for Particle seperation: %.3f cm/s\n', 100*optimal_u0);

%% Plot trajectories

% Get the screen size (screen resolution)
screenSize = get(0, 'ScreenSize');

% Define the size of each figure
figWidth = 0.5 * screenSize(3);  % 50% of screen width for the left and right sections
figHeight = (0.5 * screenSize(4)) - 60;  % 50% of screen height for each figure

% Calculate vertical spacing between figures
verticalSpacing = 50;  % Adjust as needed for space between figures

figure(1);
plot(y_grain_euler(:,1), y_grain_euler(:,2), 'b--', 'LineWidth', 3 , 'DisplayName', 'Grain Particle Trajectory(Euler)');
hold on;
plot(y_grain_rk4(:,1), y_grain_rk4(:,2), 'g-', 'LineWidth', 3 , 'DisplayName', 'Grain Particle Trajectory(RK4)');
plot(y_chaff_euler(:,1), y_chaff_euler(:,2), 'k--' , 'LineWidth', 3 , 'DisplayName', 'Chaff Particle Trajectory(Euler)');
plot(y_chaff_rk4(:,1), y_chaff_rk4(:,2), 'm-', 'LineWidth', 3 , 'DisplayName', 'Chaff Particle Trajectory(RK4)');

% Plot bin boundaries
x_bin1 = 0.1; % x-coordinate of the left boundary of Bin 1
x_bin2 = .8; % x-coordinate of the right boundary of Bin 2
x_binsep = .55; % x-coordinate of the right boundary of Bin 2
y_bin = -0.5; % y-coordinate of the bin seperation
y_ground = -.6; % y-coordinate of the ground
plot([x_bin1, x_bin1,x_binsep,x_binsep,x_binsep,x_binsep,x_bin2,x_bin2], [y_bin, y_ground,y_ground,y_bin,y_bin,y_ground, y_ground,y_bin], 'y--', 'LineWidth', 2, 'DisplayName', 'Bin Boundary'); % Boundary of Bin
plot([0,0], [-.6,.6], 'k-', 'LineWidth', 3, 'DisplayName', 'Wall'); % Boundary of wall
plot([0,0], [-.05,.05], 'r-', 'LineWidth', 3, 'DisplayName', 'Jet Inlet'); % Boundary of slot

text(.25, -.55, 'Bin: 1 ', 'FontSize', 12, 'Color', 'b');
text(.65, -.55, 'Bin: 2 ', 'FontSize', 12, 'Color', 'b');

% Plot particle positions
plot(x0, y0, 'x', 'MarkerSize', 10,'MarkerFaceColor', [0.8500 0.5250 0.0980],'MarkerEdgeColor', [0.8500 0.5250 0.0980], 'DisplayName', 'Particle Initial Position:');
plot(y_grain_rk4(end,1), y_grain_rk4(end,2), 'd', 'MarkerSize', 8,'MarkerFaceColor', [0.8500 0.5250 0.0980],'MarkerEdgeColor', [0.8500 0.5250 0.0980], 'DisplayName', 'Grain');
plot(y_chaff_rk4(end,1), y_chaff_rk4(end,2), 'p', 'MarkerSize', 8,'MarkerFaceColor', [0.8500 0.5250 0.0980],'MarkerEdgeColor', [0.8500 0.5250 0.0980], 'DisplayName', 'Chaff');

% Define the grid for x and y
x = linspace(0.0001, .75, 15); % x should start from 1 to avoid division by zero
y = linspace(-0.3, .3, 50); % y range
[X, Y] = meshgrid(x, y);

% Calculate velocity components
U = 6.2 * u0 * sqrt(h ./ X) .* exp(-50 * (Y.^2) ./ (X.^2)); % x-component of velocity
V = zeros(size(U)); % y-component of velocity (given as zero in the formula)
quiver(X, Y, U, V, 1 , 'r','DisplayName', 'Jet Velocity'); % Plot velocity vectors

% Set axis limits
xlim([-0.02 .81]);
ylim([-0.6 0.6]);

legend('Location', 'best');
xlabel('X (m)');
ylabel('Y (m)');
title('Particle Trajectories');
grid on;

set(gcf, 'Position', [figWidth, 50, figWidth, (figHeight*2)-10]);
hold off;

% Plot errors vs. dt
figure(2);
loglog(dt_values, error_euler_grain, 'o-', 'LineWidth', 3 , 'DisplayName', 'Euler Method Grain');
hold on;
loglog(dt_values, error_rk4_grain, 's-', 'LineWidth', 3 , 'DisplayName', 'RK4 Method Grain');
loglog(dt_values, error_euler_chaff, 'o-', 'LineWidth', 3 , 'DisplayName', 'Euler Method chaff');
loglog(dt_values, error_rk4_chaff, 's-', 'LineWidth', 3 , 'DisplayName', 'RK4 Method chaff');
loglog(dt_values, dt_values, 'r--', 'LineWidth', 1 , 'DisplayName', 'Ref order 1 slope');
loglog(dt_values, dt_values.^2, 'g--', 'LineWidth', 1 , 'DisplayName', 'Ref order 2 slope');
loglog(dt_values, dt_values.^3, 'b--', 'LineWidth', 1 , 'DisplayName', 'Ref order 3 slope');
loglog(dt_values, dt_values.^4, 'k--', 'LineWidth', 1 , 'DisplayName', 'Ref order 4 slope');
xlabel('Time step (dt)');
ylabel('Error');
legend('Location', 'best');
title('Error Analysis for Euler and RK4 Methods');
grid on;

set(gcf, 'Position', [figWidth/2, figHeight/2, figWidth, figHeight]);

%Plot minimum velocity for seperation

figure(3);
plot(y_chaff_optimal_rk4(:,1), y_chaff_optimal_rk4(:,2), 'c-', 'LineWidth', 4 , 'DisplayName', 'Chaff Particle (RK4)');
hold on;
% Add text to the figure
x_text_position = .1; % X coordinate for the text
y_text_position = -.3; % Y coordinate for the text
text(x_text_position, y_text_position, ['Minimum Jet Velocity for Partical Seperation: ', num2str(100*optimal_u0,4), ' cm/s'], 'FontSize', 12, 'Color', 'm');
legend show;
xlabel('X (m)');
ylabel('Y (m)');
title('optimal Jet Velocity for Particle Separation');
grid on;


% Plot bin boundaries
plot([x_bin1, x_bin1,x_binsep,x_binsep,x_binsep,x_binsep,x_bin2,x_bin2], [y_bin, y_ground,y_ground,y_bin,y_bin,y_ground, y_ground,y_bin], 'y--', 'LineWidth', 2, 'DisplayName', 'Bin Boundary'); % Boundary of Bin
plot([0,0], [-.6,.6], 'k-', 'LineWidth', 3, 'DisplayName', 'Wall'); % Boundary of wall
plot([0,0], [-.05,.05], 'r-', 'LineWidth', 3, 'DisplayName', 'Jet Inlet'); % Boundary of slot

text(.25, -.55, 'Bin: 1 ', 'FontSize', 12, 'Color', 'b');
text(.65, -.55, 'Bin: 2 ', 'FontSize', 12, 'Color', 'b');

% Plot particle positions
plot(x0, y0, 'x', 'MarkerSize', 10,'MarkerFaceColor', [0.8500 0.5250 0.0980],'MarkerEdgeColor', [0.8500 0.5250 0.0980], 'DisplayName', 'Particle Initial Position:');
plot(y_chaff_optimal_rk4(end,1), y_chaff_optimal_rk4(end,2), 'p', 'MarkerSize', 8,'MarkerFaceColor', [0.8500 0.5250 0.0980],'MarkerEdgeColor', [0.8500 0.5250 0.0980], 'DisplayName', 'Chaff');

% Define the grid for x and y
x = linspace(0.0001, .75, 15); % x should start from 1 to avoid division by zero
y = linspace(-0.3, .3, 50); % y range
[X, Y] = meshgrid(x, y);

% Calculate velocity components
U = 6.2 * u0 * sqrt(h ./ X) .* exp(-50 * (Y.^2) ./ (X.^2)); % x-component of velocity
V = zeros(size(U)); % y-component of velocity (given as zero in the formula)
quiver(X, Y, U, V, 1 , 'r','DisplayName', 'Jet Velocity'); % Plot velocity vectors

% Set axis limits
xlim([-0.02 .81]);
ylim([-0.6 0.6]);
set(gcf, 'Position', [0, 50, figWidth, (figHeight*2)-10]);
hold off;

figure(2);

%% Function to calculate drag force
function F_drag = drag_force(velocity_particle_vector, velocity_fluid_vector, diameter, rho_fluid, fluid_mu)
    Re = rho_fluid * norm(velocity_fluid_vector - velocity_particle_vector) * diameter / fluid_mu; % Reynolds number
    if Re < 800
        Cd = (24 / Re) * (1 + 0.15 * Re^0.687);
    else
        Cd = 0.44;
    end
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

%% Equations of motion
function dydt = equations_of_motion(y, mass, diameter, rho_fluid, fluid_mu, g, h, u0)
    x = y(1);
    y_pos = y(2);
    vx = y(3);
    vy = y(4);

    % Jet velocity
    u_jet = jet_velocity(x, y_pos, h, u0);

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
    fprintf('dydt', dydt);
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

%% optimal velocity

function [optimal_u0] = mim_jet_velocity(mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0, x0, y0, vx0, vy0, xc ,yc , t_span, dt, Bin_ground)
    
    h_wait = waitbar(0, 'Calculating Optimal jet Velocity for Particle seperation...');

    % Initial parameters
    u0_initial = u0; % Initial jet velocity (m/s)
    u0_min = 0; % Minimum velocity to test
    tolerance = .001; % Tolerance for convergence
    
    % Bisection method to find optimal u0
    u0_low = u0_min;
    u0_high = u0_initial;
    
    Bin_ground = -.5;
  
    % Total iterations based on tolerance
    max_iterations = ceil(log2((u0_high - u0_low) / tolerance));

    for iter = 1:max_iterations
        u0_current = (u0_low + u0_high) / 2;
        
        % Simulate light particle
        [~, y_chaff_optimal] = rk4(@(t, y) equations_of_motion(y, mass_chaff, diameter_chaff, rho_fluid, fluid_mu, g, h, u0_current), [x0, y0, vx0, vy0], t_span, dt, Bin_ground);
        
        % Check if particle falls in Bin 2
        if check_bin_landing(y_chaff_optimal, xc, yc)
            u0_high = u0_current;
        else
            u0_low = u0_current;
        end
        close(h_wait);
        progress = iter / max_iterations;
        h_wait = waitbar(progress, 'Calculating Optimal jet Velocity for Particle seperation...');
       
    end
    close(h_wait);
    optimal_u0 = u0_high;
end

%% Check if particle lands in specified bin
function lands_in_bin = check_bin_landing(trajectory, xc, yc)
    last_x = trajectory(end, 1);
    last_y = trajectory(end, 2);
    lands_in_bin = (last_x >= xc && last_y >= yc);
end

