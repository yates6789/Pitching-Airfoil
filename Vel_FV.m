% Vel_FV.m
% James Yates - 03/11/2019
% This script evaluates lamb oseen vortices

%% Variables --------------------------------------------------------------
% z - position of vortex centre in complex coordinates - [1xn]
% Circ - Circulation of vortex - [1xn]
% lo - Core radius of vortex - [1,n]
% Z - point of interest in complex coordinates = [ixj]

%% Outputs ----------------------------------------------------------------
% V - velocity at positons Z in complex conjugte form - [ixj]

function V = Vel_FV(z,Circ,lo,Z)

%% Reshaping Points of Interest -------------------------------------------
Z_reshape = reshape(Z,[numel(Z),1]);

%% Evaluating Lamb Oseen Vortex -------------------------------------------
V = sum((-1i/(2*pi))*Circ.*(1./(Z_reshape - z)).*(1 - exp(-((magnitude(Z_reshape - z)).^2)./(lo.^2))),2); 

%% Reshaping Velocity at Points of Interest -------------------------------
V = reshape(V,size(Z));

%% Removing NaN Values ----------------------------------------------------
V(isnan(V)) = 0;

end
