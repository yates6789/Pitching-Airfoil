% Vel_VS.m
% James Yates - 03/11/2019
% This script implements the Birkhoff Rott Equation for a vortex sheet

%% Variables --------------------------------------------------------------
% z - position of vortex sheet in complex coordinates - [1xn]
% Str - Strength distribution of vortex sheet - [1xn]
% Z - point of interest in complex coordinates = [ixj]

%% Outputs ----------------------------------------------------------------
% V - velocity at positons Z in complex conjugte form - [ixj]

function V = Vel_VS(z,Str,Z)

%% Reshaping Points of Interest -------------------------------------------
Z_reshape = reshape(Z,[numel(Z),1]);

%% Surface Properties -----------------------------------------------------
% Relative Vector Components
dx = real(z(2:end) - z(1:end-1));
dy = imag(z(2:end) - z(1:end-1));
% Inclination Angle
theta = atan2(dy,dx);
% Tangent Vector Components
Sx = cos(theta);
Sy = sin(theta);

%% Number of Vortex Panels ------------------------------------------------
n = 1:(length(z)-1);

%% Evaluating Berkhoff Rott Equation --------------------------------------
A = (-1i/(2*pi))*magnitude(z(n+1) - z(n))./(z(n+1) - z(n));
B = ((Z_reshape - z(n+1))./(z(n) - z(n+1))).*log((Z_reshape - z(n))./(Z_reshape - z(n+1))) + 1;
C = ((Z_reshape - z(n))./(z(n+1) - z(n))).*log((Z_reshape - z(n+1))./(Z_reshape - z(n))) + 1;
V = sum(A.*(Str(n).*B - Str(n+1).*C),2);

%% Cauchy Term ------------------------------------------------------------
% Locating Surface Points
Cauchy_Index = find(ismember(Z,z));
% Evaluating Cauchy Term
for i = 1:length(Cauchy_Index)
    Index = z == Z(Cauchy_Index(i));
    if Z(Cauchy_Index(i)) == z(end)
       V(Cauchy_Index(i)) = -0.5*Str(1)*(Sx(1) - 1i*Sy(1)) + -0.5*Str(end)*(Sx(end) - 1i*Sy(end));
    else
       V(Cauchy_Index(i)) = -0.5*Str(Index)*(Sx(Index) - 1i*Sy(Index));
    end
end

%% Integrating ------------------------------------------------------------
V = reshape(V,size(Z));

end