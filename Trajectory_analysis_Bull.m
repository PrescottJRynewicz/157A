%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz
% 504288967
% MAE 157A Trajectory Analysis Code
% Team SpaceY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; %clear all; clc;
clc
%% Input Parameters: Big Bertha

Oz_kg               = 0.0283;                               % [kg/oz]
in_m                = 0.0254;                               % [m/in]
M_dry               = 2.0*Oz_kg;                            % [kg]
M_PL                = 0.37*Oz_kg;                           % [kg]
OD                  = 1.33*in_m;                            % [m]
PD                  = 12*in_m;                              % [m]

g                   = 9.81;                                 % [m/s]
Re                  = 3.67e6;                               % []
Thrust_curve        = dlmread('B6_4.pol'); 
t_simple            = Thrust_curve(:,1);                    % [t]
t_total             = [];
T_simple            = Thrust_curve(:,2);                    %[N]
T_total             = [];


% This Loop lionearly expands the simple data points given in the G40-7 Thrust
% Profile found at http://www.thrustcurve.org/simfilesearch.jsp?id=1969 in
% order to create a better approximation of the thrust profile. ie More
% Data Points

for i = 1:length(t_simple)-1
    expanded_array_t    = linspace(t_simple(i),t_simple(i+1),100);
    expanded_array_T    = linspace(T_simple(i),T_simple(i+1),100);
    t_total             = cat(2,t_total,expanded_array_t);
    T_total             = cat(2,T_total,expanded_array_T);
end
                                  
prop_mass_init      = 5.6/1000;                             % [kg]
total_mass_motor    = 19.1/1000;                            % [kg]
struct_mass_motor   = 9.7/1000;                             % [kg]
mass_diff           = (total_mass_motor-struct_mass_motor...
                      - prop_mass_init);
M                   = total_mass_motor+M_dry+M_PL;          % [kg]
I                   = 4.90;                                 % [N-s]
Isp                 = I/(prop_mass_init*g);                                % [sec]
m_dot               = T_total./g./Isp;                      % potential correction factor of -.0004;
                                                            % this wouls be to make
                                                            % prop_mass_profile equal zero at
                                                            % end of burn.
prop_mass_profile   = [prop_mass_init];                     % [kg]
total_mass          = linspace(M,M,length(t_total));

add = 0; 
for index = 2:length(T_total)
    prop_mass_profile(index)    = prop_mass_profile(index-1) - m_dot(index)*(t_total(index)-t_total(index-1));
    total_mass(index)           = total_mass(index-1) - m_dot(index)*(t_total(index)-t_total(index-1));
    add = add + m_dot(index)*(t_total(index)-t_total(index-1));
end

% We now have Thrust, time, and mass profiles for burn time of flight. Now
% we need to calculate initial gravity, density and Drag values to begin
% loop iterations for Burn Calcs.


burn_length         = length(t_total);                      % [vec_length]
a                   = linspace(0,0,burn_length);
u                   = linspace(0,0,burn_length);
h                   = linspace(0,0,burn_length);
Af                  = pi*(OD/2)^2;                          % [m^s]
theta               = deg2rad(10);                           % [rad]
ge                  = 9.81;                                 % [m/s]
h                   = linspace(0,0,burn_length);            % [m]
g                   = linspace(0,0,burn_length);            % [m/s^2]
g(1)                = 9.81;                                 % [m/s]
rho                 = linspace(0,0,burn_length);            % [kg/m^3]
rho(1)              = 1.225;                                % [kg/m^3]
CD                  = 0.5;                                  % Initial Estimate from Stine
D                   = linspace(0,0,burn_length);            % [N]
D(1)                = 0.5*rho(1)*u(1)^2*CD*Af;              % [N]
dt_avg              = 0; 

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
    if t_total(index) >= 3.67
        index = index+1; 
        break;
    end
end 

total_mass(index)   = total_mass(index-1) - mass_diff; 
CD                  = 1.5;                                       % Will be determined from drop tests with parachute. 
Af                  = pi*(PD/2)^2;                              % [m^2]r
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

%figure; plot(t_total,u);    title('Velocity vs. Time: C_D = 0.8/2.0', 'Fontsize', 16); xlabel('Time [s]'); ylabel('Velocity [m/s]'); 
figure; plot(t_total,h);    title('Altitude vs. Time: C_D = 0.5/1.5', 'Fontsize', 16); xlabel('Time [s]'); ylabel('Altitude [m]');
% figure; plot(t_total,D);    title('Drag vs. Time', 'Fontsize', 16);     xlabel('Time [s]'); ylabel('Drag [N]');
% figure; plot(t_total,g);    title('Gravity vs. Time', 'Fontsize', 16);  xlabel('Time [s]'); ylabel('Gravity [m/s]');
% figure; plot(t_total,rho);  title('Density vs. Time', 'Fontsize', 16);  xlabel('Time [s]'); ylabel('Density [kg/m^3]');
