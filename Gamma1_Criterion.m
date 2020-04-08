% Gamma1_Criterion.m
% James Yates - 03/11/2019
% This script implements the Gamma to identify vortices

%% Variables --------------------------------------------------------------
% Z - Complex coordinates of domain - [ixj]
% V - Flow field of domain in complex form - [1xj]
% N - Number of spacing for window (must be even and > 2) - [1]
% low_lim - Lower value to define vortex core (approx. 2/pi)
% tol - Relative distance threshold between vortex centres - [1]
% Setting PLOTS = 1 will produce plots

%% Outputs ----------------------------------------------------------------
% Gamma1_cond - Gamma2 matrix with values > low_lim set to NaN - [ixj]
% Gamma1 - Gamma2 matrix - [ixj]
% Z_centre - position of vortex centres - [1xk]
% Z_core - positions assigned to each vortex centre - [(i*j)xk]
% V_core - velocity of positions asigned to each vortex core - [(i*j)xk]

function [Gamma1_cond,Gamma1,Z_centre,Z_core,V_core] = Gamma1_Criterion(Z,V,N,low_lim,tol,PLOTS)

%% Verify Number of spacing within viewing window
if mod(N,2) == 1
   msgbox('Number of spacings within viewing window must be even and > 2')
   return
end

%% Cartesian Coordinates --------------------------------------------------
X = real(Z);
Y = imag(Z);

%% Including halo ring ----------------------------------------------------
% Complex Coordinates
Z_nan = nan(size(Z) + N);
Z_nan((1+N/2):(end-N/2),(1+N/2):(end-N/2)) = Z;
% Complex Velocity
V_nan = nan(size(V) + N); 
V_nan((1+N/2):(end-N/2),(1+N/2):(end-N/2)) = V;
% Initialising Gamma2
Gamma1 = nan(size(Z) + N);

%% Gamma2 -----------------------------------------------------------------
for i = (1 + N/2):(size(Gamma1,1) - N/2)
    for j = (1 + N/2):(size(Gamma1,2) - N/2)
        % Surrounding Positions
        Z_s = Z_nan((i - N/2):(i + N/2),(j - N/2):(j + N/2)); 
        Z_s(N/2+1,N/2+1) = nan;
        % Relative Position Vectors
        dZ_s = Z_s - Z_nan(i,j); 
        dZ_s(N/2+1,N/2+1) = nan;
        % Velocity Vectors
        V_s = V_nan((i - N/2):(i + N/2),(j - N/2):(j + N/2)); 
        V_s(N/2+1,N/2+1) = nan;
        % Alpha Angles
        Alpha_s = atan2(imag(dZ_s),real(dZ_s));
        % Gamma2 - Beta Angles
        Beta_s_2 = atan2(-imag(V_s),real(V_s));    
        % Gamma2 - Theta Angles
        Theta_s_2 = Beta_s_2 - Alpha_s;    
        % Gamma2 - Sin of Angle
        Sin_s_2 = real(sin(Theta_s_2));
        % Gamma2 - Values
        Gamma1(i,j) = abs(nansum(nansum(Sin_s_2))/sum(sum(~isnan(Sin_s_2))));     
    end
end

%% Removing halo ring -----------------------------------------------------
Gamma1 = Gamma1((1+N/2):(end-N/2),(1+N/2):(end-N/2));

%% Gamma1 - Conditioning --------------------------------------------------
Gamma1_cond = Gamma1; 
Gamma1_cond(abs(Gamma1)<low_lim) = nan;

%% Gamma1 - Initial Guesses for Root Finding ------------------------------
% Cooridnates within core
x_guess = sort(X(~isnan(Gamma1_cond)));
y_guess = sort(Y(~isnan(Gamma1_cond))); 
% Removing duplicates
guess = uniquetol([x_guess,y_guess],tol,'ByRows',true);

%% Gamma1 - Root Finding --------------------------------------------------
% Function 
f = @(t) interp2(X,Y,-abs(Gamma1),t(1),t(2));
% Boundaries
lb = [min(min(X)),min(min(Y))];
ub = [max(max(X)),max(max(Y))];
% Options
options = optimoptions('fmincon','Display','off');
% Initialising Coordinate Vectors
x_roots = zeros(size(guess,1),1);
y_roots = zeros(size(guess,1),1);
for i = 1:size(guess,1)
    % Initial Guesses
    x0 = guess(i,1); 
    y0 = guess(i,2);
    % Minimiser
    Root = fmincon(f,[x0,y0],[],[],[],[],lb,ub,[],options);
    % Storing Vortex Centres
    x_roots(i) = Root(1);
    y_roots(i) = Root(2);   
end
% Removing duplicates
Root = uniquetol([x_roots,y_roots],tol,'ByRows',true);
% Complex Coordinates
Z_centre = Root(:,1)' + 1i*Root(:,2)';

%% Gamma1 - Allocating Vortex Core to Vortex Centre -----------------------
% Vortex Core
V_core = repmat(V(~isnan(Gamma1_cond)),1,length(Z_centre));                
Z_core = repmat(Z(~isnan(Gamma1_cond)),1,length(Z_centre));                 
dZ_core = Z_core - Z_centre; 
% Least Distance 
distance = sqrt((real(dZ_core).*real(dZ_core)) + (imag(dZ_core).*imag(dZ_core)));                                            
least_distance = distance == min(distance,[],2);
% Assigning core position to centre
V_core(least_distance == 0) = nan;
Z_core(least_distance == 0) = nan;

%% Gamma1 - Criterion Plots -----------------------------------------------
if PLOTS == 1
   figure(41)
   hold on, axis square
   contourf(X,Y,Gamma1)
   plot3(real(Z_centre),imag(Z_centre),ones(size(Z_centre)),'ro','MarkerFaceColor','r')
   shading interp
   title('\Gamma_1 - Criterion')
   xlabel('x - [m]')
   ylabel('y - [m]')
   xlim([min(min(X)),max(max(X))])
   ylim([min(min(Y)),max(max(Y))])
   set(gca, 'YDir','reverse')
   c = colorbar;
   title(c,'\Gamma_1')
end

end