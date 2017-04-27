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



