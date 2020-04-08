% Gamma2_Criterion.m8
% James Yates - 03/11/2019
% This script implements the Gamma to identify vortices

%% Variables --------------------------------------------------------------
% Z - Complex coordinates of domain - [ixj]
% V - Flow field of domain in complex form - [1xj]
% N - Number of spacings for window (must be even and > 2) - [1]
% low_lim - Lower value to define vortex core (approx. 2/pi)
% z_centre - first estimate of vortex centres - [1xn]
% tol - Relative distance threshold between vortex centres - [1]
% Setting PLOTS = 1 will produce plots

%% Outputs ----------------------------------------------------------------
% Gamma2_cond - Gamma2 matrix with values > low_lim set to NaN - [ixj]
% Gamma2 - Gamma2 matrix - [ixj]
% Z_centre - position of vortex centres - [1xk]
% Z_core - positions assigned to each vortex centre - [(i*j)xk]
% V_core - velocity of positions asigned to each vortex core - [(i*j)xk]
% Gamma2_core - Gamma2 of positions asigned to each vortex core - [(i*j)xk]

function [Gamma2_cond,Gamma2,Z_centre,Z_core,V_core,Gamma2_core] = Gamma2_Criterion(Z,V,N,low_lim,Z_centre,tol,PLOTS)

%% Verify Number of spacing within viewing window
if mod(N,2) == 1
   msgbox('Number of spacings within viewing window must be even and > 2')
   return
end

%% Cartesian Coordinates --------------------------------------------------
X = real(Z); Y = imag(Z);

%% Including halo ring ----------------------------------------------------
% Complex Coordinates
Z_nan = nan(size(Z) + N);
Z_nan((1+N/2):(end-N/2),(1+N/2):(end-N/2)) = Z;
% Complex Velocity
V_nan = nan(size(V) + N); 
V_nan((1+N/2):(end-N/2),(1+N/2):(end-N/2)) = V;
% Initialising Gamma2
Gamma2 = nan(size(Z) + N);

%% Gamma2 -----------------------------------------------------------------
for i = (1 + N/2):(size(Gamma2,1) - N/2)
    for j = (1 + N/2):(size(Gamma2,2) - N/2)
        % Surrounding Positions
        Z_s = Z_nan((i - N/2):(i + N/2),(j - N/2):(j + N/2)); 
        Z_s(N/2+1,N/2+1) = nan;
        % Relative Position Vectors
        dZ_s = Z_s - Z_nan(i,j); 
        dZ_s(N/2+1,N/2+1) = nan;
        % Velocity Vectors
        V_s = V_nan((i - N/2):(i + N/2),(j - N/2):(j + N/2)); 
        V_s(N/2+1,N/2+1) = nan;
        % Advective Velocity Vector
        V_adv = nansum(nansum(V_s))/sum(sum(~isnan(V_s)));
        % Corrected Velocity Vectors
        V_s = V_s - V_adv;
        % Alpha Angles
        Alpha_s = atan2(imag(dZ_s),real(dZ_s));
        % Gamma2 - Beta Angles
        Beta_s_2 = atan2(-imag(V_s),real(V_s));    
        % Gamma2 - Theta Angles
        Theta_s_2 = Beta_s_2 - Alpha_s;    
        % Gamma2 - Sin of Angle
        Sin_s_2 = real(sin(Theta_s_2));
        % Gamma2 - Values
        Gamma2(i,j) = abs(nansum(nansum(Sin_s_2))/sum(sum(~isnan(Sin_s_2))));     
    end
end

%% Removing halo ring -----------------------------------------------------
Gamma2 = Gamma2((1+N/2):(end-N/2),(1+N/2):(end-N/2));

%% Gamma2 - Conditioning --------------------------------------------------
Gamma2_cond = Gamma2;
Gamma2_cond(abs(Gamma2)<low_lim) = nan;

%% Additional Initial Guesses for Root Finding ----------------------------
% Contours
Contour = figure('visible','off');
M = contour(X,Y,Gamma2,linspace(low_lim,0.8,10));
close(Contour)
% Defining number of contours
counter = 1;
Contour_Number = 0;
while counter < size(M,2)
      counter = counter + M(2,counter) + 1;
      Contour_Number = Contour_Number + 1;
end
% Stepping through contours
counter = 1;
Z_centre_additional = zeros(1,Contour_Number);
for i = 1:Contour_Number
   % Coordinates
   x = M(1,counter+1:counter+M(2,counter)); 
   y = M(2,counter+1:counter+M(2,counter));
   % Additional Guesses
   Z_centre_additional(i) = mean(x)+1i*mean(y);
   % Increasing Counter
   counter = counter + M(2,counter) + 1;
end
Z_centre = [Z_centre,Z_centre_additional];

%% Root Finding Vortex Centre - From Gamma1 Results -----------------------
% Function 
f = @(t) interp2(X,Y,-abs(Gamma2),t(1),t(2));
% Boundaries
lb = [min(min(X)),min(min(Y))];
ub = [max(max(X)),max(max(Y))];
% Options
options = optimoptions('fmincon','Display','off');
% Initialising Coordinate Vectors
x_roots = nan(size(Z_centre));
y_roots = nan(size(Z_centre));
for i = 1:length(Z_centre)
    if ~isnan(Z_centre(i))
    % Initial Guesses
    x0 = real(Z_centre(i)); 
    y0 = imag(Z_centre(i));
    % Minimiser
    Root = fmincon(f,[x0,y0],[],[],[],[],lb,ub,[],options);
    % Storing Vortex Centres
    x_roots(i) = Root(1);
    y_roots(i) = Root(2);   
    end
end
% Removing duplicates
Root = uniquetol([x_roots',y_roots'],tol,'ByRows',true);
% Complex Coordinates
Z_centre = Root(:,1)' + 1i*Root(:,2)';
% Removing Nan
Z_centre(isnan(Z_centre)) = [];

%% Removing Points without Centre -----------------------------------------
% Contours
Contour = figure('visible','off');
M = contour(X,Y,Gamma2,[low_lim,low_lim]);
close(Contour)
% Defining number of contours
counter = 1;
Contour_Number = 0;
while counter < size(M,2)
      counter = counter + M(2,counter) + 1;
      Contour_Number = Contour_Number + 1;
end
% Stepping through contours
counter = 1;
for i = 1:Contour_Number
   % Coordinates
   x = M(1,counter+1:counter+M(2,counter))'; 
   y = M(2,counter+1:counter+M(2,counter))';  
   % Removing Points without centre
   if sum(inpolygon(real(Z_centre),imag(Z_centre),x,y)) == 0
      Gamma2_cond(inpolygon(real(Z),imag(Z),x,y)) = nan;
   end
   counter = counter + M(2,counter) + 1;
end

%% Gamma2 - Allocating Vortex Core to Centre - Assuming no axisymmetry ----
% Vortex Core
V_core = repmat(V(~isnan(Gamma2_cond)),1,length(Z_centre));                
Z_core = repmat(Z(~isnan(Gamma2_cond)),1,length(Z_centre));   
Gamma2_core = repmat(Gamma2(~isnan(Gamma2_cond)),1,length(Z_centre));
dZ_core = Z_core - Z_centre; 
% Least Distance 
distance = sqrt((real(dZ_core).*real(dZ_core)) + (imag(dZ_core).*imag(dZ_core)));                                       
least_distance = distance == min(distance,[],2);
% Assigning core position to centre
V_core(least_distance == 0) = nan;
Z_core(least_distance == 0) = nan;
Gamma2_core(least_distance == 0) = nan;
% Gamma2 averaged centre
Z_centre = nanmean(Z_core.*Gamma2_core,1)./nanmean(Gamma2_core,1);

%% Gamma2 - Centre Values for Plots ---------------------------------------
% Evaluating Gamma2 Values at vortex centres
Gamma2_centre = zeros(size(Z_centre));
for i = 1:length(Z_centre)
    Gamma2_centre(i) = interp2(X,Y,Gamma2,real(Z_centre(i)),imag(Z_centre(i)));
end

%% Gamma2 - Criterion Plots -----------------------------------------------
if PLOTS == 1
   figure(51)
   hold on, axis square
   contourf(X,Y,Gamma2)
   plot3(real(Z_centre),imag(Z_centre),Gamma2_centre,'ro','MarkerFaceColor','r')
   shading interp
   title('\Gamma_2 - Criterion')
   xlabel('x - [m]')
   ylabel('y - [m]')
   xlim([min(min(X)),max(max(X))])
   ylim([min(min(Y)),max(max(Y))])
   set(gca, 'YDir','reverse')
   c = colorbar;
   title(c,'\Gamma_2')
   caxis([0,1])
end

end