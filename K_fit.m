% K_fit.m
% James Yates - 04/03/2020
% This script minimises flow field residuals to relax the kutta condition 

%% Variables --------------------------------------------------------------
% z - Complex coordinates of vortex sheet (must be closed loop) - [1xn]
% Vk - known velocities at airfoil surface in complex form - [1xn]
% U - free stream Velocity - [1]
% V_without_sheet - velocity without vortex sheet in complex form - [ixj]
% Z - Complex coordinates of domain - [ixj]
% V - Desired flow field for comparison - [ixj]

%% Outputs ----------------------------------------------------------------
% K - optimal net-strength at the trailing edge - [1]
% R_Residual - R squared when applying the optimal value - [1]
% R_Kutta - R squared when applying the Kutta condition - [1]

function [K,R_Residual,R_Kutta] = K_fit(z,Vk,U,V_without_sheet,Z,V)

%% Tolerances for finite difference approximation of optima ---------------
h = 1e-4;
tol = 1e-6;

%% Vortex Panel Method - Function -----------------------------------------
Str_f = @(t) Str_VS(z,Vk,U,t,"NO");

%% Velocity Field - Function ----------------------------------------------
V_f = @(t) Vel_VS(z,Str_f(t),Z) + V_without_sheet;

%% Sum of Squares of Residuals - Function ---------------------------------
SS_res = @(t) nansum(nansum((magnitude(V_f(t)) - magnitude(V)).^2));

%% Total Sum of Squares - Function ----------------------------------------
SS_tot = nansum(nansum(magnitude(V).^2));

%% R Squared - Function ---------------------------------------------------
R_f = @(t) 1 - SS_res(t)/SS_tot;

%% R Squared First Derivative - Function ----------------------------------
dR_f = @(t) (R_f(t + h) - R_f(t))/h;

%% R Squared Second Derivative - Function ---------------------------------
d2R_f = @(t) (R_f(t + h) - 2*R_f(t) + R_f(t - h))/(h*h);

%% Least Square Minimisation ----------------------------------------------
% Initial Guess as Kutta Condition
K = 0;
% Initial Error
ea = 1;
% Secant Method
while ea > tol
    % Previous Guess
    K_0 = K;
    % First Derivative
    dR = dR_f(K);
    % Second Derivative
    d2R = d2R_f(K);
    % Next Guess
    K = K - dR/d2R;
    % Approximate Error
    ea = abs((K-K_0)/K);
end

%% R Squared Value --------------------------------------------------------
% Residual Minimisation
R_Residual = R_f(K);
% Kutta Condition
R_Kutta = R_f(0);

end