%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz
% 504288967
% MAE 157A Trajectory Analysis Code
% Team SpaceY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Constants
Oz_kg               = 0.0283;                               % [kg/oz]
in_m                = 0.0254;                               % [m/in]
g                   = 9.81;                                 % [m/s]
Re                  = 3.67e6;                               % []

%% Total Vehicle Input Parameters
total_mass          = 0.7;                                  % [kg]
CD_coast_recovery   = [0.75 1.5];                           % [ ]
original_D          = 1.33*in_m;                            % [m]
recovery_D          = 18*in_m;                              % [m]
launch_angle        = 5;                                    % [deg]

%% Motor Inputs 
prop_mass           = 37.9/1000;                            % [kg]
total_mass_motor    = 83.9/1000;                            % [kg]
final_motor_mass    = 38.4/1000;                            % [kg]
total_impulse       = 76.83;                                % [N-s]
deploy_time         = 6.05; 
Thrust_curve        = dlmread('Thrust_curve.pol');          % [N]

%% Interpolate Thrust curve for high accuracy in numerical integration. 
t_simple            = Thrust_curve(:,1);                    % [t]
t_total             = [];                                   % [t]    
T_simple            = Thrust_curve(:,2);                    % [N]
T_total             = [];                                   % [N]
for i = 1:length(t_simple)-1
    expanded_array_t    = linspace(t_simple(i),t_simple(i+1),100);
    expanded_array_T    = linspace(T_simple(i),T_simple(i+1),100);
    t_total             = cat(2,t_total,expanded_array_t);
    T_total             = cat(2,T_total,expanded_array_T);
end

%% Vehicle Calcs
final_mass          = total_mass-total_mass_motor ...       % [ ] 
                      + final_motor_mass;                   % [kg]
Isp                 = total_impulse/(prop_mass*g);          % [sec]
m_dot               = T_total./g./Isp;                      % [kg/s]
prop_mass_profile   = linspace(0,0,length(t_total));        % [kg]
prop_mass_profile(1)=prop_mass;                             % [kg]
total_mass          = linspace(total_mass,total_mass,...    % [kg]
                      length(t_total));                     % [kg]
                  
%% mass profile from changing prop mass
for index = 2:length(T_total)
    prop_mass_profile(index)    = prop_mass_profile(index-1) - m_dot(index)*(t_total(index)-t_total(index-1));
    total_mass(index)           = total_mass(index-1) - m_dot(index)*(t_total(index)-t_total(index-1));
end


%% Setup all variables for numerical calcs. 
burn_length         = length(t_total);                      % [vec_length]
a                   = linspace(0,0,burn_length);            % [m/s^2]
u                   = linspace(0,0,burn_length);            % [m/s]
h                   = linspace(0,0,burn_length);            % [m]
Af                  = pi*(original_D/2)^2;                  % [m^s]
theta               = deg2rad(launch_angle);                % [rad]
ge                  = 9.81;                                 % [m/s]
h                   = linspace(0,0,burn_length);            % [m]
g                   = linspace(0,0,burn_length);            % [m/s^2]
g(1)                = 9.81;                                 % [m/s]
rho                 = linspace(0,0,burn_length);            % [kg/m^3]
rho(1)              = 1.225;                                % [kg/m^3]
CD                  = 0.8;                                  % Initial Estimate from Stine
D                   = linspace(0,0,burn_length);            % [N]
D(1)                = 0.5*rho(1)*u(1)^2*CD*Af;              % [N]
dt_avg              = 0; 


%% Loop 1: Calculate all values from T-0 to end of motor burn. 
% This could also be integrated into the mass profile loop, 
% but code already very quickly, so not currently a concern. 
for index = 1:(burn_length-1)
    dt              = (t_total(index+1)-t_total(index));
    dt_avg          = dt_avg + dt; 
    a(index+1)      = (T_total(index) - D(index))/total_mass(index) - g(index)*cos(theta);
    u(index+1)      = u(index) + (T_total(index)/total_mass(index) - D(index)/total_mass(index) ...
                                - g(index)*cos(theta))*dt;
    h(index+1)      = h(index) + u(index)*dt*cos(theta);
    g(index+1)      = ge*(Re/(Re+h(index+1)));
    rho(index+1)    = real(1.2*exp(-2.9*10^-5*h(index+1)^1.15));
    D(index+1)      = 0.5*rho(index+1)*u(index+1)^2*CD*Af; 
end


%% Loop 2: Coast time from end of burn to parachute deploy
%  This loop assumes that the parachute deploys before rocket 
% reaches max apogee due to pure coast. 
dt                  = dt_avg/burn_length;                   % [s] New Time step from average thrust dt.
for index = burn_length:10000000
    total_mass(index)= total_mass(index-1);
    
    a(index+1)      = -D(index)/total_mass(index) - g(index)*cos(theta);
    u(index+1)      = u(index) + (-D(index)/total_mass(index)*dt ...
                                - g(index)*cos(theta)*dt);
    h(index+1)      = h(index) + u(index)*dt*cos(theta);
    g(index+1)      = ge*(Re/(Re+h(index+1)));
    rho(index+1)    = real(1.2*exp(-2.9*10^-5*h(index+1)^1.15));
    D(index+1)      = 0.5*rho(index+1)*u(index+1)^2*CD*Af; 
    t_total(index+1)= t_total(index)+dt; 
    if t_total(index) >= deploy_time
        index = index+1; 
        break;
    end
end 


%% Loop 3: Update values to include parachute. 
%  runs until velocity = 0 and vehicle begins to fall
% back to earth. 
total_mass(index)   = final_mass; 
CD                  = 2.0;                                       % Will be determined from drop tests with parachute. 
Af                  = pi*(recovery_D/2)^2;                              % [m^2]r
D(index)            = 0.5*rho(index)*u(index)^2*CD*Af; 

for index = index:100000000
    total_mass(index)= total_mass(index-1);
    
    u(index+1)      = u(index) + (-D(index)/total_mass(index)*dt ...
                                - g(index)*dt);
    h(index+1)      = h(index) + u(index)*dt*cos(theta);
    g(index+1)      = ge*(Re/(Re+h(index+1)));
    rho(index+1)    = real(1.2*exp(-2.9*10^-5*h(index+1)^1.15));
    D(index+1)      = 0.5*rho(index+1)*u(index+1)^2*CD*Af; 
    t_total(index+1)= t_total(index)+dt; 
    if u(index+1) <= 0 
        break;
    end
end


%% Loop 4: Just switch direction of drag force
for index = index:100000000
    total_mass(index)= total_mass(index-1);
    
    u(index+1)      = u(index) + (D(index)/total_mass(index)*dt ...
                                - g(index)*dt);
    h(index+1)      = h(index) + u(index)*dt*cos(theta);
    g(index+1)      = ge*(Re/(Re+h(index+1)));
    rho(index+1)    = real(1.2*exp(-2.9*10^-5*h(index+1)^1.15));
    D(index+1)      = 0.5*rho(index+1)*u(index+1)^2*CD*Af; 
    t_total(index+1)= t_total(index)+dt; 
    if h(index+1) <= 0 
        break;
    end
end

%% Plots

figure; plot(t_total,u);    title('Velocity vs. Time', 'Fontsize', 16); xlabel('Time [s]'); ylabel('Velocity [m/s]'); 
figure; plot(t_total,h);    title('Altitude vs. Time', 'Fontsize', 16); xlabel('Time [s]'); ylabel('Altitude [m]');
% figure; plot(t_total,D);    title('Drag vs. Time', 'Fontsize', 16);     xlabel('Time [s]'); ylabel('Drag [N]');
% figure; plot(t_total,g);    title('Gravity vs. Time', 'Fontsize', 16);  xlabel('Time [s]'); ylabel('Gravity [m/s]');
% figure; plot(t_total,rho);  title('Density vs. Time', 'Fontsize', 16);  xlabel('Time [s]'); ylabel('Density [kg/m^3]');
