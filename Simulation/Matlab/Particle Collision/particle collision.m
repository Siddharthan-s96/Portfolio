%% Programme: Collisions of particles, comparison of hard-sphere and DEM 
%
% Author: Siddharthan Somasundaram ( matriculation number : 254304 )
%
% Brief description: This program models the collision dynamics of two particles using both hard-sphere 
% and soft-sphere (Discrete Element Method or DEM) approaches.It investigates the differences in collision 
% modeling by comparing an idealized, instantaneous collision (hard-sphere) with a more realistic, continuous
% interaction (soft-sphere). The program utilizes analytical solutions for the hard-sphere case to determine
% post-collision velocities based on the coefficient of restitution. For the soft-sphere approach, it employs
% a numerical integration scheme (Leapfrog method) to simulate the particles' motion, incorporating a linear 
% spring model to represent the contact force during collisions. By comparing the results obtained from both 
% methods, including contact points, post-collision velocities, and collision duration, the program 
% provides insights into the limitations and advantages of each approach in modeling particle interactions.
%
% Submission date: 03/02/2025
%
%%

% Clear workspace and command window
clear variables
close all
clc

% Display header for program start
fprintf('%s\n', repmat('*', 1, 140)); 
fprintf('%s', repmat('*', 1, 54)); fprintf(' PARTICLE COLLISION SIMULATION '); fprintf('%s\n', repmat('*', 1, 55));
fprintf('%s', repmat('*', 1, 58)); fprintf(' TIME INTEGRATION START '); fprintf('%s\n', repmat('*', 1, 58));
fprintf('%s\n\n', repmat('*', 1, 140)); 

%% Simulation Parameters

% Particle 1 initial conditions
p1_init_pos = [0, 0, 0];             % [x, y, z] coordinates
p1_init_vel = [0, 0, 0];             % [vx, vy, vz] components
p1_diameter = 1.0;                    % Units: meters
p1_mass = 0.05;                       % Units: kg

% Particle 2 initial conditions
p2_init_pos = [1.1, 1.3, 0];         % [x, y, z] coordinates
p2_init_vel = [-1.0, -1.0, 0];       % [vx, vy, vz] components
p2_diameter = 1.0;                    % Units: meters
p2_mass = 0.05;                       % Units: kg

% Collision parameters
coeff_restitution_elastic = 1.0;      % Perfectly elastic collision
coeff_restitution_inelastic = 0.75;   % Partially elastic collision

% Time integration parameters
timestep_vis = 0.01;                  % Larger timestep for visualization (seconds)
timestep_dem = 1e-4;                  % Smaller timestep for DEM simulation (seconds)
time_total = 2;                       % Total simulation time (seconds)

%% Collision Detection and Analysis

t1 = tic; % Start the first timer

% Calculate time and positions of collision based on initial conditions
[time_collision] = calculate_collision_params( ...
    p1_init_pos, p1_init_vel, p2_init_pos, p2_init_vel, p1_diameter, p2_diameter);

[p1_pre_collision, p2_pre_collision] = calculate_position_euler(p1_init_pos, p1_init_vel, p2_init_pos, p2_init_vel, time_collision);

if time_collision > 0
    % Calculate initial contact point at the moment of collision
    contact_point_init = (p1_pre_collision + p2_pre_collision) / 2;

    % Compute post-collision velocities for elastic collision
    disp('Computing elastic collision using Hard-sphere algorithm...');
    [p1_vel_hs_elastic, p2_vel_hs_elastic] = compute_hardsphere_velocity( ...
        p1_pre_collision, p1_init_vel, p2_pre_collision, p2_init_vel, ...
        p1_mass, p2_mass, coeff_restitution_elastic);

    % Compute post-collision velocities for inelastic collision
    disp('Computing inelastic collision using Hard-sphere algorithm...');
    [p1_vel_hs_inelastic, p2_vel_hs_inelastic] = compute_hardsphere_velocity( ...
        p1_pre_collision, p1_init_vel, p2_pre_collision, p2_init_vel, ...
        p1_mass, p2_mass, coeff_restitution_inelastic);

    % Set collision duration to zero for instantaneous hard-sphere model
    duration_hs = 0;
    
    time_post_collision = time_collision + duration_hs;
    
    [p1_pos_collision, p2_pos_collision] = calculate_position_euler(p1_init_pos, p1_init_vel, p2_init_pos, p2_init_vel, time_post_collision);


    contact_point_final = (p1_pos_collision + p2_pos_collision) / 2;

else
    % Display message if no collision occurs within simulation timeframe
    disp('No collision detected in simulation timeframe.');
end

elapsed1 = toc(t1);

%% Soft-sphere DEM Simulation

t2 = tic;

% Compute elastic collision using the Soft-sphere DEM algorithm
disp('Computing elastic collision using Soft-sphere algorithm...');
[p1_pos_ss_elastic, p1_vel_ss_elastic, p2_pos_ss_elastic, p2_vel_ss_elastic, ...
 overlap_ss_elastic, normal_ss_elastic] = compute_ss_leapfrog(p1_init_pos, p1_init_vel, ...
    p2_init_pos, p2_init_vel, p1_mass, p2_mass, p1_diameter, p2_diameter, ...
    timestep_dem, time_total, coeff_restitution_elastic);

% Compute inelastic collision using the Soft-sphere DEM algorithm
disp('Computing inelastic collision using Soft-sphere algorithm...');
[p1_pos_ss_inelastic, p1_vel_ss_inelastic, p2_pos_ss_inelastic, p2_vel_ss_inelastic, ...
 overlap_ss_inelastic, normal_ss_inelastic] = compute_ss_leapfrog(p1_init_pos, p1_init_vel, ...
    p2_init_pos, p2_init_vel, p1_mass, p2_mass, p1_diameter, p2_diameter, ...
    timestep_dem, time_total, coeff_restitution_inelastic);



%% Analysis of Collision Results

% Find the contact point at the start of the collision (soft-sphere model)
overlap_data = overlap_ss_inelastic;                            % Overlap data for the inelastic case
idx_contact_start = find(overlap_data > 0, 1, 'first');         % Index of the first contact
p1_pos_contact = p1_pos_ss_elastic(idx_contact_start, :);       % Position of particle 1 at the start of collision
p2_pos_contact = p2_pos_ss_elastic(idx_contact_start, :);       % Position of particle 2 at the start of collision
contact_point_ss_init = (p1_pos_contact + p2_pos_contact) / 2;  % Initial contact point (average of positions)

% Extract final velocities for elastic collisions (hard-sphere and soft-sphere)
p1_vel_final_hs_elastic = p1_vel_hs_elastic;                 % Final velocity of particle 1 (hard-sphere, elastic)
p2_vel_final_hs_elastic = p2_vel_hs_elastic;                 % Final velocity of particle 2 (hard-sphere, elastic)
p1_vel_final_ss_elastic = p1_vel_ss_elastic(end, :);         % Final velocity of particle 1 (soft-sphere, elastic)
p2_vel_final_ss_elastic = p2_vel_ss_elastic(end, :);         % Final velocity of particle 2 (soft-sphere, elastic)

% Extract final velocities for inelastic collisions (hard-sphere and soft-sphere)
p1_vel_final_hs_inelastic = p1_vel_hs_inelastic;             % Final velocity of particle 1 (hard-sphere, inelastic)
p2_vel_final_hs_inelastic = p2_vel_hs_inelastic;             % Final velocity of particle 2 (hard-sphere, inelastic)
p1_vel_final_ss_inelastic = p1_vel_ss_inelastic(end, :);     % Final velocity of particle 1 (soft-sphere, inelastic)
p2_vel_final_ss_inelastic = p2_vel_ss_inelastic(end, :);     % Final velocity of particle 2 (soft-sphere, inelastic)

% Calculate the duration of the collision for the soft-sphere model
idx_collision_start = find(overlap_ss_inelastic > 0, 1, 'first');           % Start index of collision
idx_collision_end = find(overlap_ss_inelastic > 0, 1, 'last');              % End index of collision
duration_ss = (idx_collision_end - idx_collision_start) * timestep_dem;     % Duration of collision (time)

% Calculate the final contact point for inelastic collision (soft-sphere model)
p1_pos_final = p1_pos_ss_inelastic(idx_collision_end, :);    % Final position of particle 1
p2_pos_final = p2_pos_ss_inelastic(idx_collision_end, :);    % Final position of particle 2
contact_point_ss_final = (p1_pos_final + p2_pos_final) / 2;  % Final contact point (average of positions)

elapsed2 = toc(t2);

%% Results Display

% Table formatting parameters
col_width = [7, 60, 15, 25, 25]; % Column widths
create_line = @(char, width) repmat(char, 1, width); % Function to create a separator line
center_text = @(text, width) sprintf('%*s%s%*s', ...
    floor((width - length(text)) / 2), '', text, ceil((width - length(text)) / 2), ''); % Center text
format_str = sprintf('|%%%ds|%%%ds|%%%ds|%%%ds|%%%ds|\n', ...
    col_width(1), col_width(2), col_width(3), col_width(4), col_width(5)); % Table format string

% Display results header
fprintf('\n%s\n', repmat('#', 1, 140));
fprintf('\nResults:\n');

% Display collision time
time_text = sprintf('%.2f Seconds', time_collision);
fprintf('\n%s %s %50s\n\n', center_text('1', col_width(1)), 'Collision occurrence time', center_text(time_text, col_width(3)));

% Print table headers and separators
fprintf(format_str, ...
    create_line('-', col_width(1)), create_line('-', col_width(2)), create_line('-', col_width(3)), ...
    create_line('-', col_width(4)), create_line('-', col_width(5)));

fprintf(format_str, ...
    center_text('S.No', col_width(1)), center_text('Parameter', col_width(2)), center_text('Particle', col_width(3)), ...
    center_text('Hard-Sphere Model', col_width(4)), center_text('Soft-Sphere Model', col_width(5)));

fprintf(format_str, ...
    create_line('-', col_width(1)), create_line('-', col_width(2)), create_line('-', col_width(3)), ...
    create_line('-', col_width(4)), create_line('-', col_width(5)));

% Row a: Initial contact points
fprintf(format_str,center_text('a)', col_width(1)), center_text('Initial Contact Point', col_width(2)), ...
    center_text('-', col_width(3)),center_text(sprintf('[%.3f, %.3f, %.3f]', contact_point_init), ...
    col_width(4)),center_text(sprintf('[%.3f, %.3f, %.3f]', contact_point_ss_init), col_width(5)));

% Print separator
fprintf(format_str,create_line('-', col_width(1)), create_line('-', col_width(2)), ...
    create_line('-', col_width(3)),create_line('-', col_width(4)), create_line('-', col_width(5)));

% Row b: Elastic collision results
fprintf(format_str, '', '', center_text('Particle 1', col_width(3)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p1_vel_final_hs_elastic), col_width(4)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p1_vel_final_ss_elastic), col_width(5)));

fprintf(format_str,center_text('b)', col_width(1)), ...
    center_text('Elastic Post-collision velocities (CoR = 1.0)', col_width(2)), '', '', '');

fprintf(format_str, '', '', center_text('Particle 2', col_width(3)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p2_vel_final_hs_elastic), col_width(4)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p2_vel_final_ss_elastic), col_width(5)));

% Print separator
fprintf(format_str, ...
    create_line('-', col_width(1)), create_line('-', col_width(2)), create_line('-', col_width(3)), ...
    create_line('-', col_width(4)), create_line('-', col_width(5)));

% Row c: Inelastic collision results
fprintf(format_str, '', '', center_text('Particle 1', col_width(3)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p1_vel_final_hs_inelastic), col_width(4)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p1_vel_final_ss_inelastic), col_width(5)));

fprintf(format_str, center_text('c)', col_width(1)) ...
    , center_text('Inelastic Post-collision velocities (CoR = 0.75)', col_width(2)), '', '', '');

fprintf(format_str, '', '', center_text('Particle 2', col_width(3)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p2_vel_final_hs_inelastic), col_width(4)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', p2_vel_final_ss_inelastic), col_width(5)));

% Print separator
fprintf(format_str, ...
    create_line('-', col_width(1)), create_line('-', col_width(2)), create_line('-', col_width(3)), ...
    create_line('-', col_width(4)), create_line('-', col_width(5)));

% Row d: Collision duration
duration_hs_text = sprintf('%.4f Seconds', duration_hs);
duration_ss_text = sprintf('%.4f Seconds', duration_ss);
fprintf(format_str,center_text('d)', col_width(1)), center_text('Duration of collision', ...
    col_width(2)), center_text('-', col_width(3)), center_text(duration_hs_text, col_width(4)), ...
    center_text(duration_ss_text, col_width(5)));

% Print separator
fprintf(format_str, ...
    create_line('-', col_width(1)), create_line('-', col_width(2)), create_line('-', col_width(3)), ...
    create_line('-', col_width(4)), create_line('-', col_width(5)));

% Row e: Final contact points
fprintf(format_str, ...
    center_text('e)', col_width(1)), center_text('Final Contact Point', col_width(2)), center_text('-', col_width(3)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', contact_point_final), col_width(4)), ...
    center_text(sprintf('[%.3f, %.3f, %.3f]', contact_point_ss_final), col_width(5)));

% Print final separator
fprintf(format_str, ...
    create_line('-', col_width(1)), create_line('-', col_width(2)), create_line('-', col_width(3)), ...
    create_line('-', col_width(4)), create_line('-', col_width(5)));

fprintf('\n%s\n', repmat('#', 1, 140));


%% Visualization

fprintf('\nComputational time for Hard sphere model: %.4f seconds\n', elapsed1);
fprintf('Computational time for Soft sphere model: %.4f seconds\n\n', elapsed2);

% Calculate the factor (how much slower/faster the soft-sphere model is)
if elapsed1 > 0  % Avoid division by zero
    factor = elapsed2 / elapsed1;
    fprintf('Soft-sphere model is %.2f times slower than the hard-sphere model.\n\n', factor);
elseif elapsed2 > 0
    factor = elapsed1 / elapsed2;
        fprintf('Hard-sphere model is %.2f times slower than the soft-sphere model.\n\n', factor);
else
    fprintf('Both models had zero computational time. Cannot calculate factor.\n\n');
end


% Visualize Hard Sphere - Elastic Case
disp('Visualizing elastic collision (Hard-sphere model)...');
visualize_collision_hs(p1_init_pos, p1_init_vel, p2_init_pos, p2_init_vel, ...
    p1_diameter, p2_diameter, p1_vel_hs_elastic, p2_vel_hs_elastic, ...
    timestep_vis, time_total, time_collision, coeff_restitution_elastic, duration_hs);

% Visualize Hard Sphere - Inelastic Case
disp('Visualizing inelastic collision (Hard-sphere model)...');
visualize_collision_hs(p1_init_pos, p1_init_vel, p2_init_pos, p2_init_vel, ...
    p1_diameter, p2_diameter, p1_vel_hs_inelastic, p2_vel_hs_inelastic, ...
    timestep_vis, time_total, time_collision, coeff_restitution_inelastic, duration_hs);

% Visualize Soft Sphere - Elastic Case
disp('Visualizing elastic collision (Soft-sphere model)...');
visualize_collision_ss(p1_pos_ss_elastic, p1_vel_ss_elastic, p2_pos_ss_elastic, p2_vel_ss_elastic, ...
    p1_diameter, p2_diameter, timestep_dem, coeff_restitution_elastic, idx_contact_start, duration_ss);

% Visualize Soft Sphere - Inelastic Case
disp('Visualizing inelastic collision (Soft-sphere model)...');
visualize_collision_ss(p1_pos_ss_inelastic, p1_vel_ss_inelastic, p2_pos_ss_inelastic, p2_vel_ss_inelastic, ...
    p1_diameter, p2_diameter, timestep_dem, coeff_restitution_inelastic, idx_contact_start, duration_ss);

%% Function to compute post-collision velocities using the Hard-Sphere model

% This function calculates the velocities of two colliding particles after
% the collision, based on their masses, initial positions, and velocities.
% It supports both elastic (e = 1.0) and inelastic (0 < e < 1) collisions.

function [v1_post, v2_post] = compute_hardsphere_velocity(pos1, vel1, pos2, vel2, m1, m2, e)
    
        % Compute post-collision velocities for inelastic collision using impulse-based method
        
        % Calculate the normal unit vector between the two particles
        n = (pos1 - pos2) / norm(pos1 - pos2);
        
        % Calculate the relative velocity of the two particles
        v12 = vel1 - vel2;

        % Step 3: Compute the dot product
        v12_dot_n = dot(v12, n);
        
        % Step 4: Compute the reduced mass
        m_star = (m1 * m2) / (m1 + m2);
        
        % Step 5: Compute the impulse
        j = -(1 + e) * v12_dot_n * m_star * n;

        % Update post-collision velocities based on the impulse
        v1_post = vel1 + (j / m1) ; % Post-collision velocity of particle 1
        v2_post = vel2 - (j / m2) ; % Post-collision velocity of particle 2
   
end

%% Function to calculate collision parameters (time and positions)
% This function determines the time of collision (`t_col`) and the positions
% of two particles (`pos1_col` and `pos2_col`) at the moment of collision.
% Inputs:
% - pos1, vel1: Initial position and velocity of particle 1
% - pos2, vel2: Initial position and velocity of particle 2
% - diameter1, diameter2: Diameters of particle 1 and particle 2
% Outputs:
% - t_col: Time of collision (-1 if no collision occurs)
% - pos1_col, pos2_col: Positions of particle 1 and particle 2 at collision

function [t_col] = calculate_collision_params(pos1, vel1, pos2, vel2, diameter1, diameter2)
    % Compute relative position and velocity vectors
    x12 = pos1 - pos2; % Relative position vector
    v12 = vel1 - vel2; % Relative velocity vector
    
    % Compute collision distance (sum of radii of the two particles)
    collision_distance = 0.5 * (diameter1 + diameter2);
    
    % Compute dot products for collision detection
    x12_dot_v12 = dot(x12, v12); % Dot product of relative position and velocity
    v12_dot_v12 = dot(v12, v12); % Dot product of relative velocity with itself
    x12_dot_x12 = dot(x12, x12); % Dot product of relative position with itself
    
    % Check if a valid collision occurs by evaluating the discriminant
    sqrt_term = (x12_dot_v12)^2 - v12_dot_v12 * (x12_dot_x12 - collision_distance^2);
    if sqrt_term < 0 || v12_dot_v12 == 0
        % If discriminant is negative or particles have no relative motion, no collision occurs
        t_col = -1; % No collision
        pos1_col = []; % Empty positions
        pos2_col = [];
        return;
    end
    
    % Compute the two possible collision times
    t1 = (-x12_dot_v12 - sqrt(sqrt_term)) / v12_dot_v12; % First possible collision time
    t2 = (-x12_dot_v12 + sqrt(sqrt_term)) / v12_dot_v12; % Second possible collision time
    
    % Select the smallest positive collision time
    if t1 >= 0 && t2 >= 0
        t_col = min(t1, t2); % Both times are positive, take the earlier one
    elseif t1 >= 0
        t_col = t1; % Only t1 is valid
    elseif t2 >= 0
        t_col = t2; % Only t2 is valid
    else
        % If neither time is positive, no collision occurs
        t_col = -1;
        
        return;
    end
    
end

%%

function [pos1_col, pos2_col] = calculate_position_euler(pos1, vel1, pos2, vel2, t_col)
      
    if t_col <= 0 
        pos1_col = [];
        pos2_col = [];
   
    else
        % Compute the positions of both particles at the time of collision
        pos1_col = pos1 + vel1 * t_col; % Position of particle 1 at collision
        pos2_col = pos2 + vel2 * t_col; % Position of particle 2 at collision        
    end

end

%%

% Performs a leapfrog integration to simulate the collision of two particles.
function [positions1, velocities1, positions2, velocities2, overlaps, normals] = compute_ss_leapfrog(position1, velocity1, position2, velocity2, ...
                                                                               mass1, mass2, diameter1, diameter2, ...
                                                                               time_step, total_time, CoR)
    % Number of time steps
    num_steps = ceil(total_time / time_step);
    
    % Preallocate arrays for positions and velocities
    positions1 = zeros(num_steps, 3); % Particle 1 positions
    velocities1 = zeros(num_steps, 3); % Particle 1 velocities
    positions2 = zeros(num_steps, 3); % Particle 2 positions
    velocities2 = zeros(num_steps, 3); % Particle 2 velocities
    % Preallocate arrays for overlap and normal
    overlaps = zeros(num_steps, 1);    % Overlap at each timestep
    normals = zeros(num_steps, 3);     % Normal vector at each timestep
    % Initialize positions and velocities
    positions1(1, :) = position1;
    velocities1(1, :) = velocity1;
    positions2(1, :) = position2;
    velocities2(1, :) = velocity2;
       
    % Loop over time steps
    for i = 1:num_steps-1
        % Compute overlap and normal vector
        [overlaps(i), normals(i,:)] = compute_overlap_and_normal(positions1(i, :), positions2(i, :), diameter1, diameter2);
        if overlaps(i) > 0
            % Compute forces during collision
            [force1, force2] = compute_collision_forces(velocities1(i, :), velocities2(i, :), ...
                                                         CoR, overlaps(i), normals(i,:));
        else
            % No collision, no forces
            force1 = [0, 0, 0];
            force2 = [0, 0, 0];
        end
        % Step 1: Compute initial acceleration
        acceleration1 = force1 / mass1;
        acceleration2 = force2 / mass2;
        
        % Step 2: Kick (Predictor) - Half step velocity update
        velocity1_half = velocities1(i, :) + 0.5 * acceleration1 * time_step;
        velocity2_half = velocities2(i, :) + 0.5 * acceleration2 * time_step;
        
        % Step 3: Drift - Update positions
        positions1(i+1, :) = positions1(i, :) + velocity1_half * time_step;
        positions2(i+1, :) = positions2(i, :) + velocity2_half * time_step;
        
        [overlaps(i+1), normals(i+1,:)] = compute_overlap_and_normal(positions1(i+1, :), positions2(i+1, :), diameter1, diameter2);
        if overlaps(i+1) > 0
            % Compute forces during collision
            [force1, force2] = compute_collision_forces(velocity1_half, velocity2_half, ...
                                                         CoR, overlaps(i+1), normals(i+1,:));
        else
            % No collision, no forces
            force1 = [0, 0, 0];
            force2 = [0, 0, 0];
        end
        % Step 5: Compute new acceleration
        acceleration1 = force1 / mass1;
        acceleration2 = force2 / mass2;
        
        % Step 6: Kick (Corrector) - Half step velocity update
        velocities1(i+1, :) = velocity1_half + 0.5 * acceleration1 * time_step;
        velocities2(i+1, :) = velocity2_half + 0.5 * acceleration2 * time_step;
    end
end

%%

function [force1, force2] = compute_collision_forces(vel1, vel2, CoR, overlap, normal)
  
    % Spring constant for loading (controls material hardness)
    k_loading = 750; % Adjust as needed for material properties
    
    % Compute relative velocity along the normal direction
    relative_velocity = vel2 - vel1; % Velocity of particle 2 relative to particle 1
    vn = dot(relative_velocity, normal); % Normal component of relative velocity
    
    % Determine if loading or unloading
    if vn < 0 % Loading phase (particles are approaching each other)
        spring_constant = k_loading;
    else % Unloading phase (particles are separating)
        k_unloading = CoR^2 * k_loading; % Scale spring constant using CoR
        spring_constant = k_unloading;
    end
    
    % Compute normal force
    repulsive_force_magnitude = max(0, spring_constant * overlap); % Ensure no attractive force
    
    % Compute force vectors
    force1 = -repulsive_force_magnitude * normal; % Force on particle 1
    force2 = +repulsive_force_magnitude * normal; % Force on particle 2
end

%% COMPUTE_OVERLAP_AND_NORMAL Computes the overlap distance and normal vector between two circles (2D) or spheres (3D).

function [overlap, normal] = compute_overlap_and_normal(pos1, pos2, diameter1, diameter2)
    
    % Compute the distance vector between the centers of the two particles
    distance_vector = pos2 - pos1;
    
    % Compute the distance (magnitude of the distance vector)
    distance = norm(distance_vector);
    
    % Compute the collision distance (sum of radii)
    collision_distance = 0.5 * (diameter1 + diameter2);
    
    % Compute overlap (scalar value)
    overlap = max(0, collision_distance - distance);
    
    % Compute the normal vector (unit vector pointing from particle 1 to particle 2)
    tolerance = 1e-6; % Small tolerance for floating-point comparisons
    if distance > tolerance
        normal = distance_vector / distance;
    else
        % If the distance is very small (near perfect overlap), set normal to zero
        normal = [0; 0; 0]; % Corrected to column vector for consistency
    end
end

%%

function visualize_collision_hs(pos1, vel1, pos2, vel2, diameter1, diameter2, ...
                             vel1_post, vel2_post, dt, t_limit, t_col,cor,duration_hs)

    % Time vector
    t = 0:dt:t_limit;

    % Preallocate particle trajectories
    traj1 = zeros(length(t), 3);
    traj2 = zeros(length(t), 3);

    % Compute trajectories
    for i = 1:length(t)
        if t(i) <= t_col
            traj1(i, :) = pos1 + vel1 * t(i);
            traj2(i, :) = pos2 + vel2 * t(i);
        else
            traj1(i, :) = pos1 + vel1 * t_col + vel1_post * (t(i) - t_col);
            traj2(i, :) = pos2 + vel2 * t_col + vel2_post * (t(i) - t_col);
        end
    end

    % Get the screen size (screen resolution)
    screenSize = get(0, 'ScreenSize');

    % Define the size of each figure
    figWidth = 0.5 * screenSize(3);  % 50% of screen width
    figHeight = (0.5 * screenSize(4)) - 60;  % 50% of screen height

    % Set up figure
    figure;
    hold on;
    grid on;
    legend('show');

    if cor == 1
        set(gcf, 'Position', [0, 50, figWidth, ((figHeight * 2) - 10)]); 
        title('Hard-Sphere elastic Particle Collision Visualization');
    else
        set(gcf, 'Position', [figWidth, 50, figWidth, ((figHeight * 2) - 10)]);
        title('Hard-Sphere in-elastic Particle Collision Visualization');
    end
   
    % Fixed axis limits
    xlim([-5, 5]);
    ylim([-5, 5]);

    xlabel('X position (m)');
    ylabel('Y position (m)');
    
    % Create sphere visualization
    [x_sphere, y_sphere, z_sphere] = sphere(20);
    
    % Define RGB colors 
    particle1_color = [72, 209, 204] / 255;   % Convert to MATLAB scaled RGB
    particle2_color = [240, 128, 128] / 255; % Convert to MATLAB scaled RGB
    
    % Define Particle 1 with Teal color
    particle1 = surf(zeros(2), zeros(2), zeros(2), ...
        'FaceColor', particle1_color, 'EdgeColor', 'none', 'DisplayName', 'Particle 1');
    
    % Define Particle 2 with Coral color
    particle2 = surf(zeros(2), zeros(2), zeros(2), ...
        'FaceColor', particle2_color, 'EdgeColor', 'none', 'DisplayName', 'Particle 2');
    
    

    % Add dynamic text display for time and velocities
    time_text = text(0, -3, '', 'FontSize', 12, 'VerticalAlignment', 'top');

    % Create text objects for table with proper alignment
    collision_time_text = text(-4.5, 4.5, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                                'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    
    collision_Duration_text = text(-4.5, 4.25, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                                'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    % Table header with proper spacing
    header_text = sprintf('%-15s %-25s %-25s', '', 'Pre-Collision Velocity', 'Post-Collision Velocity');

    table_header = text(-4.5, 3.75, header_text, 'FontSize', 10, 'FontName', 'Courier New', ...
                        'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    
    % Horizontal line
    line_text = repmat('-', 1, 65);
    table_line = text(-4.5, 3.5, line_text, 'FontSize', 10, 'FontName', 'Courier New', ...
                      'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    
    % Table rows with proper spacing
    table_row1 = text(-4.5, 3.25, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                      'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    table_row2 = text(-4.5, 3, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                      'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');

    % Animation loop
    for i = 1:length(t)
        % Update particle positions
        set(particle1, 'XData', 0.5 * diameter1 * x_sphere + traj1(i, 1), ...
                      'YData', 0.5 * diameter1 * y_sphere + traj1(i, 2), ...
                      'ZData', 0.5 * diameter1 * z_sphere + traj1(i, 3));
        set(particle2, 'XData', 0.5 * diameter2 * x_sphere + traj2(i, 1), ...
                      'YData', 0.5 * diameter2 * y_sphere + traj2(i, 2), ...
                      'ZData', 0.5 * diameter2 * z_sphere + traj2(i, 3));
        
        % Update text for current time and velocity
        if t(i) <= t_col
            % Pre-collision text
            set(time_text, 'String', sprintf('Time: %.3f s\n\nVelocity ( m / s )\n     Particle 1 : [%.3f, %.3f, %.3f]\n     Particle 2 : [%.3f, %.3f, %.3f]', ...
                t(i), vel1(1), vel1(2), vel1(3), vel2(1), vel2(2), vel2(3)));
        else
            % Update text for post-collision time and velocities
            set(time_text, 'String', sprintf('Time: %.3f s\n\nVelocity ( m / s )\n     Particle 1 : [%.3f, %.3f, %.3f]\n     Particle 2 : [%.3f, %.3f, %.3f]', ...
                t(i), vel1_post(1), vel1_post(2), vel1_post(3), vel2_post(1), vel2_post(2), vel2_post(3)));

            % Dynamically align table
                collision_time_str = sprintf('Collision Time        : %.3f   Seconds', t_col);
            collision_Duration_str = sprintf('Duration of collision : %.4f Seconds', duration_hs);
            col_width_1 = 15; % Width for first column
            col_width_2 = 25; % Width for second column
            col_width_3 = 25; % Width for third column

            header_text = sprintf('%-*s%-*s%-*s', col_width_1, '', col_width_2, 'Pre-Collision Velocity', col_width_3, 'Post-Collision Velocity');
            line_text = repmat('-', 1, col_width_1 + col_width_2 + col_width_3);

            row1_text = sprintf('%-*s%-*s%-*s', ...
                col_width_1, 'Particle 1:', col_width_2, sprintf('[%.3f, %.3f, %.3f]', vel1(1), vel1(2), vel1(3)), ...
                col_width_3, sprintf('[%.3f, %.3f, %.3f]', vel1_post(1), vel1_post(2), vel1_post(3)));

            row2_text = sprintf('%-*s%-*s%-*s', ...
                col_width_1, 'Particle 2:', col_width_2, sprintf('[%.3f, %.3f, %.3f]', vel2(1), vel2(2), vel2(3)), ...
                col_width_3, sprintf('[%.3f, %.3f, %.3f]', vel2_post(1), vel2_post(2), vel2_post(3)));

            % Display table
            set(collision_time_text, 'String', collision_time_str, 'Visible', 'on');
            set(collision_Duration_text, 'String', collision_Duration_str, 'Visible', 'on');
            set(table_header, 'String', header_text, 'Visible', 'on');
            set(table_line, 'String', line_text, 'Visible', 'on');
            set(table_row1, 'String', row1_text, 'Visible', 'on');
            set(table_row2, 'String', row2_text, 'Visible', 'on');
       
        end

        drawnow;
        pause(0.01);
    end
    grid off;
    hold off;
end

%%

function visualize_collision_ss(positions1, velocities1, positions2, velocities2, ...
                                  diameter1, diameter2, dt, cor,contact_start, duration_ss)

    % Time vector
    num_steps = size(positions1, 1);
    time = (0:num_steps) * dt;
    
    ss_t_col = dt*contact_start ;

    % Get the screen size (screen resolution)
    screenSize = get(0, 'ScreenSize');

    % Define the size of each figure
    figWidth = 0.5 * screenSize(3);  % 50% of screen width
    figHeight = (0.5 * screenSize(4)) - 60;  % 50% of screen height

    % Set up figure
    figure;
    hold on;
    grid on;
    legend('show');

    if cor == 1
        set(gcf, 'Position', [0, 50, figWidth, ((figHeight * 2) - 40)]); 
        title('Soft-Sphere elastic Particle Collision Visualization');
    else
        set(gcf, 'Position', [figWidth, 50, figWidth, ((figHeight * 2) - 40)]);
        title('Soft-Sphere in-elastic Particle Collision Visualization');
    end    

    % Fixed axis limits
    xlim([-5, 5]);
    ylim([-5, 5]);

    xlabel('X position (m)');
    ylabel('Y position (m)');


    % Create sphere visualization
    [x_sphere, y_sphere, z_sphere] = sphere(20);
    
    % Define RGB colors 
    particle1_color = [72, 209, 204] / 255;   % Convert to MATLAB scaled RGB
    particle2_color = [240, 128, 128] / 255; % Convert to MATLAB scaled RGB
    
    % Define Particle 1 with Teal color
    particle1 = surf(zeros(2), zeros(2), zeros(2), ...
        'FaceColor', particle1_color, 'EdgeColor', 'none', 'DisplayName', 'Particle 1');
    
    % Define Particle 2 with Coral color
    particle2 = surf(zeros(2), zeros(2), zeros(2), ...
        'FaceColor', particle2_color, 'EdgeColor', 'none', 'DisplayName', 'Particle 2');
    
    
    % Add dynamic text display for time and velocities
    time_text = text(0, -3, '', 'FontSize', 12, 'VerticalAlignment', 'top');

    % Create text objects for table with proper alignment
    collision_time_text = text(-4.5, 4.5, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                                'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    
    collision_Duration_text = text(-4.5, 4.25, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                                'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');

    % Table header with proper spacing
    header_text = sprintf('%-15s %-25s %-25s', '', 'Pre-Collision Velocity', 'Post-Collision Velocity');
    table_header = text(-4.5, 3.75, header_text, 'FontSize', 10, 'FontName', 'Courier New', ...
                        'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    
    % Horizontal line
    line_text = repmat('-', 1, 65);
    table_line = text(-4.5, 3.5, line_text, 'FontSize', 10, 'FontName', 'Courier New', ...
                      'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    
    % Table rows with proper spacing
    table_row1 = text(-4.5, 3.25, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                      'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');
    table_row2 = text(-4.5, 3, '', 'FontSize', 10, 'FontName', 'Courier New', ...
                      'HorizontalAlignment', 'left', 'Color', [0.0, 0.0, 0.55], 'Visible', 'off');



   % Animation loop with skipping frames
    frame_skip = 100; % Plot every 10th frame
    for i = 1:frame_skip:num_steps
        % Update particle positions
        set(particle1, 'XData', 0.5 * diameter1 * x_sphere + positions1(i, 1), ...
                       'YData', 0.5 * diameter1 * y_sphere + positions1(i, 2), ...
                       'ZData', 0.5 * diameter1 * z_sphere + positions1(i, 3));
        set(particle2, 'XData', 0.5 * diameter2 * x_sphere + positions2(i, 1), ...
                       'YData', 0.5 * diameter2 * y_sphere + positions2(i, 2), ...
                       'ZData', 0.5 * diameter2 * z_sphere + positions2(i, 3));
        
        if time(i) <= ss_t_col
             % Pre-collision text
            set(time_text, 'String', sprintf('Time: %.3f s\n\nVelocity ( m / s )\n     Particle 1 : [%.3f, %.3f, %.3f]\n     Particle 2 : [%.3f, %.3f, %.3f]', ...
                time(i), velocities1(i, 1), velocities1(i, 2), velocities1(i, 3), velocities2(i, 1), velocities2(i, 2), velocities2(i, 3)));

        else
            set(time_text, 'String', sprintf('Time: %.3f s\n\nVelocity ( m / s )\n     Particle 1 : [%.3f, %.3f, %.3f]\n     Particle 2 : [%.3f, %.3f, %.3f]', ...
                time(i), velocities1(i, 1), velocities1(i, 2), velocities1(i, 3), velocities2(i, 1), velocities2(i, 2), velocities2(i, 3)));
        
            % Dynamically align table
                collision_time_str = sprintf('Collision Time        : %.3f   Seconds', ss_t_col);
            collision_Duration_str = sprintf('Duration of collision : %.4f Seconds', duration_ss);
            col_width_1 = 15; % Width for first column
            col_width_2 = 25; % Width for second column
            col_width_3 = 25; % Width for third column

            header_text = sprintf('%-*s%-*s%-*s', col_width_1, '', col_width_2, 'Pre-Collision Velocity', col_width_3, 'Post-Collision Velocity');
            line_text = repmat('-', 1, col_width_1 + col_width_2 + col_width_3);

            row1_text = sprintf('%-*s%-*s%-*s', ...
                col_width_1, 'Particle 1:', col_width_2, sprintf('[%.3f, %.3f, %.3f]', velocities1(1, 1), velocities1(1, 2), velocities1(1, 3)), ...
                col_width_3, sprintf('[%.3f, %.3f, %.3f]', velocities1(i, 1), velocities1(i, 2), velocities1(i, 3)));

            row2_text = sprintf('%-*s%-*s%-*s', ...
                col_width_1, 'Particle 2:', col_width_2, sprintf('[%.3f, %.3f, %.3f]', velocities2(1, 1), velocities2(1, 2), velocities2(1, 3)), ...
                col_width_3, sprintf('[%.3f, %.3f, %.3f]', velocities2(i, 1), velocities2(i, 2), velocities2(i, 3)));

            % Display table
            set(collision_time_text, 'String', collision_time_str, 'Visible', 'on');
             set(collision_Duration_text, 'String', collision_Duration_str, 'Visible', 'on');
            set(table_header, 'String', header_text, 'Visible', 'on');
            set(table_line, 'String', line_text, 'Visible', 'on');
            set(table_row1, 'String', row1_text, 'Visible', 'on');
            set(table_row2, 'String', row2_text, 'Visible', 'on');

        
        end
    
        drawnow;
        pause(0.01); % Adjust if needed
    end
    grid off;
    hold off;
end

