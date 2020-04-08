% Circulation.m
% James Yates - 11/01/2020
% This script performs a closed loop intergral to determine circulation
% from the results of the Gamma2 criterion

%% Variables --------------------------------------------------------------
% Z_centre - position of vortex centres - [1]
% Z_core - positions assigned to each vortex centre - [nx1]
% V_core - velocity of positions asigned to each vortex core - [nx1]
% Gamma2_core - Gamma2 of positions asigned to each vortex core - [nx1]
% Setting PLOTS = 1 will produce plots

%% Outputs ----------------------------------------------------------------
% Circ_LO - circulation of the fit vortex - [1]
% Radius_LO - core radius of the fit vortex - [1]
% Circ_LO_Error - error in circulation (relative) - [1]
% Radius_LO_Error - error in radius (relative) - [1]
% R_squared_LO - R squared of fit - [1]

function [Circ_LO,Radius_LO,Circ_LO_Error,Radius_LO_Error,R_squared_LO] = Circulation(Z_centre,Z_core,V_core,Gamma2_core,PLOTS)

%% Centering Vortex Core --------------------------------------------------
Z_core = Z_core - Z_centre;
Z_core(isnan(V_core)) = [];
Gamma2_core(isnan(V_core)) = [];
V_core(isnan(V_core)) = [];

%% Meshing Vortex Core ----------------------------------------------------
% Domain
X = unique(sort(real(Z_core))); Y = unique(sort(imag(Z_core)));
% Meshing Domain
[X,Y] = meshgrid(X,Y); Z_core_2 = X + 1i*Y;
% Find active points of Domain
Gamma2_core_2 = nan(size(Z_core_2));
V_core_2 = nan(size(Z_core_2));
for i = 1:length(Z_core)
    Gamma2_core_2(Z_core_2 == Z_core(i)) = Gamma2_core(i);
    V_core_2(Z_core_2 == Z_core(i)) = V_core(i);
end

%% Return function if closed loop is not possible for any radius - 1st Run
if numel(unique(X)) <= 2 || numel(unique(Y)) <= 2
   % Lamb Oseen Vortex
   Circ_LO = 0;
   Radius_LO = 0;
   Circ_LO_Error = 0;
   Radius_LO_Error = 0;
   R_squared_LO = 0;
   % Return Function
   return
end

%% Circulation across Vortex Core -----------------------------------------
% Contours
Contour = figure('visible','off');
M = contour(X,Y,Gamma2_core_2,unique(Gamma2_core));
close(Contour)
% Defining number of contours
counter = 1;
Contour_Number = 0;
while counter < size(M,2)
      counter = counter + M(2,counter) + 1;
      Contour_Number = Contour_Number + 1;
end
% Initialising vectors
Gamma2 = nan(1,Contour_Number);
Circ_R = nan(1,Contour_Number);
R = nan(1, Contour_Number);
% Stepping through contours
counter = 1;   % contour level index
for i = 1:Contour_Number
   % Coordinates
   x = M(1,counter+1:counter+M(2,counter))'; 
   y = M(2,counter+1:counter+M(2,counter))';
   Z = x + 1i*y;
   if (Z(1) == Z(end)) && (numel(Z) > 3) && inpolygon(0,0,x,y)  
      % Ensuring Counter Clockwise Direction
      K = convhull(x,y);
      % Coordinates
      x = x(K);
      y = y(K);      
      % Velocity
      u = interp2(X,Y,real(V_core_2),x,y);
      v = interp2(X,Y,-imag(V_core_2),x,y);
      % Relative Vector Components
      dx = x(2:end) - x(1:end-1);
      dy = y(2:end) - y(1:end-1);
      dl = sqrt(dx.*dx + dy.*dy);
      % Inclination Angle
      theta = atan2(dy,dx);
      % Tangent Vector Components
      tnx = cos(theta);
      tny = sin(theta);
      % Tangent Velocity
      dCirc = (u(1:end-1).*tnx + v(1:end-1).*tny);
      % Circulation Intergral
      Circ_R(i) = sum(dCirc.*dl);
      % Contour Level
      Gamma2(i) = M(1,counter);
      % Radius 
      R(i) = mean(sqrt(x.^2 + y.^2));       
   end
   counter = counter + M(2,counter) + 1;
end
% Removing NaNs
R(isnan(Circ_R)) = [];
Circ_R(isnan(Circ_R)) = [];
% Sorting Values 
SORT = sortrows([R',Circ_R']);
R = SORT(:,1);
Circ_R = SORT(:,2);

%% Return function if closed loop is not possible for any radius - 2nd Run
if numel(R) < 3
   % Lamb Oseen Vortex
   Circ_LO = 0;
   Radius_LO = 0;
   Circ_LO_Error = 0;
   Radius_LO_Error = 0;
   R_squared_LO = 0;
   % Return Function
   return
end

%% Lamb Oseen Vortex - Direct Fit -----------------------------------------
% Estimating Start Point
A = mean(Circ_R);
B = max(R);
% Upper and Lower Bounds
if mean(Circ_R) > 0
   lb = [0,0];
   ub = [10*A,B];
else
   lb = [10*A,0];
   ub = [0,B];
end
% Fitting Circulation Distribution
[Fit_LO,Stats_LO] = fit(R,Circ_R,fittype('a*(1 - exp(-(x/b)^2))'),'StartPoint',[A,B],'lower',lb,'upper',ub);
% Circulation
Circ_LO = Fit_LO.a;
Radius_LO = abs(Fit_LO.b);
% Error Calculation - 95%
Confidence_LO = confint(Fit_LO);
Circ_LO_Error = abs((Confidence_LO(2,1)-Circ_LO)/Circ_LO);
Radius_LO_Error = abs((Confidence_LO(2,2)-Radius_LO)/Radius_LO);
if isnan(Radius_LO_Error)
   Radius_LO_Error = 0;
end
% R Squared
R_squared_LO = Stats_LO.rsquare;

%% Plotting Root Finding --------------------------------------------------
if PLOTS == 1
   figure(61), hold on, grid on
   % Lamb Oseen Plots
   plot(R,Circ_R,'bo')
   plot(linspace(0,max(R),100),Circ_LO*(1 - exp(-(linspace(0,max(R),100)./Radius_LO).^2)),'r-','LineWidth',2)
   xlabel('R - [m]'), xlim([0,max(R)])
   ylabel('\Gamma - [m^2/s]')
   title('Lamb Oseen Vortex')
   legend('Data','Direct Fit')
end

end
