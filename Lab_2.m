%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz
% 504288967
% MAE 157A Lab 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Constant

XF                  = 0.68;                     % [m]
CR                  = 0.055;                    % [m]
CT                  = 0.027;                    % [m]
XS                  = 0.041;                    % [m]
S                   = 0.05;                     % [m]
LN                  = 0.162;                    % [m]
D                   = 0.0338;                   % [m]
CN_alpha_n          = 2.0;                      % [1/rad]

%% Nose Cone Cp
alpha_n             = 0.500;                    % [rad]
x_n                 = alpha_n*LN;               % [m*rad]

%% Main Stage Fins

K_fb                = 1 + (D/(2*S + D));        % [ ]
LF                  = sqrt(S^2 + (XS + CT/2  ...% [ ]
                      -CR/2)^2);                % [m]
beta                = 16;                       % [ ]
CN_alpha_f          = beta*(S/D)^2/(1 + sqrt(...% [ ]
                     1 + (2*LF/(CR+CT))^2));    % [ ]
CN_alpha_fb         = K_fb*CN_alpha_f;          % [ ]
x_f                 = XF + (XS*(CR+2*CT))/   ...% [ ]
                      (3*CR+CT) + 1/6*(CR +  ...% [ ]
                      CT - CR*CT/(CR+CT));      % [m]
%% Final CP Calc

CN_alpha_T          = CN_alpha_n + CN_alpha_fb; % [ ]
x_cp                = (CN_alpha_n*x_n +     ... % [ ]
                      CN_alpha_fb*x_f)/     ... % [ ]
                      CN_alpha_T;               % [m]
x_cg                = 0.508;                    % [m]
SM                  = (x_cp - x_cg)/D*100;      % [%]



%% Team Tasks

clark_y             = dlmread('Lab2_data.lvm');
sd8020              = dlmread('Lab2_data_2.lvm'); 
rectangle           = dlmread('Lab2_data_3.lvm');


Drag                = 1;                        % drag
Lift                = 2;                        % Lift
Side                = 3;                        % side
Mx                  = 4; 
My                  = 5; 
Mz                  = 6; 
vel                 = 7; 

S                   = 0.03097;                  % [m^2]
rho                 = 1.225;                    % [kg/m^3] 
L                   = 0.2032;                   % [m]
mu                  = 1.789E-5;                 % [kg/m/s]
Re_clark            = rho.*L.*clark_y(:,vel)/mu;% [ ] 
Re_sd8020           = rho.*L.*sd8020(:,vel)/mu; % [ ] 
Re_rect             = rho.*L.*rectangle(:,vel)...
                      /mu;                      % [ ] 
                  
%% Drag Coefficient Vs. Reynolds Number

CD_clark            = 2*clark_y(:,Drag)./(rho.*clark_y(:,vel).^2.*S);
CD_sd8020           = 2*sd8020(:,Drag)./(rho.*sd8020(:,vel).^2.*S);
CD_rect             = 2*rectangle(:,Drag)./(rho.*rectangle(:,vel).^2.*S);


figure; plot(Re_clark, CD_clark); title('Drag Coefficient vs. Reynolds Number: Clark-Y Airfoil'); xlabel('Reynolds Number'); ylabel('Drag Coefficient'); 
figure; plot(Re_sd8020, CD_sd8020); title('Drag Coefficient vs. Reynolds Number: SD8020 Airfoil'); xlabel('Reynolds Number'); ylabel('Drag Coefficient'); 
figure; plot(Re_rect, CD_rect); title('Drag Coefficient vs. Reynolds Number: Rectangular Airfoil'); xlabel('Reynolds Number'); ylabel('Drag Coefficient'); 
figure; plot(Re_clark, clark_y(:,Drag)); title('Drag vs. Reynolds Number: Clark-Y Airfoil'); xlabel('Reynolds Number'); ylabel('Drag [N]'); 
figure; plot(Re_sd8020, sd8020(:,Drag)); title('Drag vs. Reynolds Number: SD8020 Airfoil'); xlabel('Reynolds Number'); ylabel('Drag [N]'); 
figure; plot(Re_rect, rectangle(:,Drag)); title('Drag vs. Reynolds Number: Rectangular Airfoil'); xlabel('Reynolds Number'); ylabel('Drag [N]'); 

%% Flutter

b                   = .05;                      % [m]
AR                  = 6^2/(8*6);                % [ ]
gamma               = 1.4;                      % [ ]
R                   = 287;                      % [J/kg/K]
T                   = 298;                      % [ ]
a                   = sqrt(gamma*R*T);          % [m/s]
G                   = 2.275E9;                  % [Pa]

%








