%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Individual Calculations                           %
% Written and Developed By:                         %
% Janelle T. Rogers, UCLA MAE Undergraduate Student %
% Date:                                             %
% 09 April 2017                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%
clc; clear all; close all;
 
%% iv.
 
sigma_yield = 1.450*10^8;                           % yield stress [N/m^2]
OD_in = 2;                                          % outer-diameter [in]
OD = 2*0.0254;                                      % outer-diameter [m] USE THIS ONE
wt_in = 0.035;                                      % wall-thickness [in]
wt = 0.035*0.0254;                                  % wall-thickness [m] USE THIS ONE
r_outer = OD/2;                                     % outer radius [m]
r_inner = (OD - (2 * wt))/2;                        % inner radius [m]
A_outer = pi*r_outer^2;                             % outer-tube area [m^2]
A_inner = pi*r_inner^2;                             % inner-tube area [m^2]
A = A_outer - A_inner;                              % cross-sectional area [m^2]
F_crityield = sigma_yield*A;                        % yield force [N]
 
n = 4;                                              % value of end when both ends are fixed
E = 6.9*(10^10);                                    % modulus of elasticity [N/m^2], Solidworks material properties value
I = (pi/4)*(r_outer^4 - r_inner^4);                 % area moment of interia of hollow cylinder [m^4]
L_in = 6;                                           % tube length [in]
L = 6*0.0254;                                       % tube length [m] USE THIS ONE
F_critbuck = (n*(pi^2)*E*I)/(L^2);                  % buckling force [N]
 
%% v.
 
Lnew = sqrt((n*(pi^2)*E*I)/(sigma_yield*A));        % length required for buckling [m]
