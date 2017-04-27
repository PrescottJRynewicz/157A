%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz
% 504288967
% MAE 157A Lab 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 


min_stress          = [4.83E7 2.49E7 5.20E7]./1000;         % [kN/m^2]
max_stress          = [1.99E8 9.65E7 5.20E7]./1000;         % [kN/m^2]
yield_strength      = [1.45E8 1.45E8 8.27E8]./1000;         % [kN/m^2]
thickness           = [0.035 0.070 0.035].*0.0254;          % [m]
E                   = [6.9E10 6.9E10 113.8E9]/1000;         % [Pa : kN/m^2]
min_displacement    = [0 0 0];                              % [m]
max_displacement    = [3.17E-1 1.61E-1 2.08E-1].*0.0254;    % [m]
min_FOS             = [0.73 1.50 4.28];                     % [ ]
rho                 = [2.8 2.8 4.43]./1000.*100^3;          % [kg/m^3]
stw                 = zeros(1,3);                           % [ ]
A                   = zeros(1,3);                           % [m^2]
F_cr_yield          = zeros(1,3);                           % [kN]
F_cr_buckle         = zeros(1,3);                           % [kN]
I                   = zeros(1,3);                           % [m^4]
L                   = [6 6 6]*0.0254;                       % [m]
D                   = [2 2 2].*.0254; 



%% Strength to wieght Ratio

for index = 1:3
    A(index)                = pi*((D(index)/2)^2 - ((D(index) - 2*thickness(index))/2)^2); 
    stw(index)              = yield_strength(index)./rho(index); 
    F_cr_yield(index)       = yield_strength(index)*A(index);
    I(index)                = pi/64*((0.0508)^4 - (0.0508 - 2*thickness(index))^4); 
    F_cr_buckle(index)      = 4*pi^2*E(index)*I(index)/L(index)^2; 
    L(index)                = sqrt(4*pi^2*E(index)*I(index)/(yield_strength(index)*A(index)));    
end


%% Team Tasks

max_force           = 1171.55;                              % [N]
applied_force       = 110;                                  % [N]
FOS                 = max_force/applied_force;              % [%]


