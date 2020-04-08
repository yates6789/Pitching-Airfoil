% Pitching_Motion_Characterisation.m
% James Yates - 23/01/2020
% This script characterises the flow field of a pitching airfoil using
% based on the analysis of experimental data
% - External vortices are identifies using Gamma1 and Gamma2 criteria and
%   detected vortices are modelled by fitting a Lamb Oseen Vortex
% - Airfoil Surface is modelled using the Vortex Panel Method
% - The first model applies the Kutta condition at the trailing edge
% - The second model compares the estimated solution to the experimental
%   data to define and optimal trailinge edge condition
clc; clear; close all;

%% Vortex Detection Settings ----------------------------------------------
Gamma_tol = 0.08;                                                           % Relative distance threshold allowed between vortex centres
Gamma_lim = 0.5;                                                            % Lower value to define vortex core (approx. 2/pi)
Gamma1_Target = 0.1;                                                        % Target detection ratio for Gamma1 criteria
Gamma2_Target = 0.64;                                                       % Target detection ratio for Gamma2 criteria
N1 = 2;                                                                     % Number of spacings for Gamma1 detection window (must be even and > 2) - [1]
N2 = 16;                                                                    % Number of spacings for Gamma2 detection window (must be even and > 2) - [1]

%% Flow Field Characteristics ---------------------------------------------
U = 0.215;                                                                  % Free Stream Velocity
Chord_Length = 0.15;                                                        % Chord Length of Airfoil
Pitching_Axis = (106.98 + 144.45*1i)/1000;                                  % Position of Pitching Axis in complex coordinates
theta = deg2rad(62.3121);                                                   % Angle of Attack - [rad]
omega = deg2rad(16.8386);                                                   % Pitching Rate - [rad/s]

%% Airfoil Contour --------------------------------------------------------
NACA = [0,0,0.18];
PANELS = 100;
z = NACA_Airfoil(NACA(1),NACA(2),NACA(3),Chord_Length,Pitching_Axis,theta,PANELS,'NO');

%% Importing Data ---------------------------------------------------------
% Raw Data
Data = readmatrix('test.csv');
% Coordinates
x = Data(:,1);
y = Data(:,2);
% Velocity
u = Data(:,3);
v = Data(:,4);

%% Reshaping vectors to arrays --------------------------------------------
% Coordinates
x = (reshape(x,[numel(unique(x)),numel(unique(y))])')/1000;
y = (reshape(y,[numel(unique(x)),numel(unique(y))])')/1000;
% Velocity
u = reshape(u,[numel(unique(x)),numel(unique(y))])';
v = reshape(v,[numel(unique(x)),numel(unique(y))])';

%% Complex Forms ----------------------------------------------------------
% Coordinates
Z = x + 1i*y;
% Velocity
V = u - 1i*v;

%% Gamma1 - Criterion -----------------------------------------------------
[Gamma1_cond,Gamma1,Z1_centre,Z1_core,V1_core] = Gamma1_Criterion(Z,V-U,N1,Gamma_lim,Gamma_tol,"NO");

%% Gamma2 - Criterion -----------------------------------------------------
[Gamma2_cond,Gamma2,Z2_centre,Z2_core,V2_core,Gamma2_core] = Gamma2_Criterion(Z,V,N2,Gamma_lim,Z1_centre,Gamma_tol,"NO");

%% Vortex Fitting ---------------------------------------------------------
% Initialising vortex fitting matrices
Circ = zeros(size(Z2_centre));
Radius = zeros(size(Z2_centre));
Circ_Error = zeros(1,length(Z2_centre));
Radius_Error = zeros(1,length(Z2_centre));
R_squared = zeros(size(Z2_centre));
% Vortex Fitting
for i = 1:length(Z2_centre)
    [Circ(i),Radius(i),Circ_Error(i),Radius_Error(i),R_squared(i)] = ...
     Circulation(Z2_centre(:,i),Z2_core(:,i),V2_core(:,i),Gamma2_core(:,i),"NO");
end

%% Removing Vortices outlier vortices -------------------------------------
Average_Error = 0.5*(Circ_Error + Radius_Error);
STD_Error = std(Average_Error(Average_Error > 0));
Mean_Error = mean(Average_Error(Average_Error > 0));
R_squared((Average_Error > (Mean_Error + 1.96*STD_Error))|(Average_Error > 1)) = 0;

%% Removing Vortices without a fit ----------------------------------------
Z2_centre(R_squared <= 0) = [];
Circ(R_squared <= 0) = [];
Radius(R_squared <= 0) = [];
Circ_Error(R_squared <= 0) = [];
Radius_Error(R_squared <= 0) = [];
R_squared(R_squared <= 0) = [];

%% Removing Vortices within Airfoil ---------------------------------------
Circ(inpolygon(real(Z2_centre),imag(Z2_centre),real(z),imag(z))) = [];
Radius(inpolygon(real(Z2_centre),imag(Z2_centre),real(z),imag(z))) = [];
Circ_Error(inpolygon(real(Z2_centre),imag(Z2_centre),real(z),imag(z))) = [];
Radius_Error(inpolygon(real(Z2_centre),imag(Z2_centre),real(z),imag(z))) = [];
R_squared(inpolygon(real(Z2_centre),imag(Z2_centre),real(z),imag(z))) = [];
Z2_centre(inpolygon(real(Z2_centre),imag(Z2_centre),real(z),imag(z))) = [];

%% Detection Ratio --------------------------------------------------------
% Viewing Window Length
L_1 = (N1/2)*(x(1,2)-x(1,1));
L_2 = (N2/2)*(x(1,2)-x(1,1));
% Detection Ratio Approximation
Ratio_1 = L_1./Radius;
Ratio_2 = L_2./Radius;

%% Next Guess for Gamma1 Window Size --------------------------------------
% Window Length Guess
L1_i = Gamma1_Target*Radius;
% Window Size Guess
N1_i = 2*ceil(2*L1_i/(x(1,2)-x(1,1))/2);

%% Next Guess for Gamma2 Window Size --------------------------------------
% Window Length Guess
L2_i = Gamma2_Target*Radius;
% Window Size Guess
N2_i = 2*ceil(2*L2_i/(x(1,2)-x(1,1))/2);

%% known velocity at midpoint of vortex panels ----------------------------
% Midpoints
z_m = [0.5*(z(2:end)+z(1:end-1)),z(1)];
% Free vortices
V_FV_k = Vel_FV(Z2_centre,Circ,Radius,z_m);
% Pitching Motion
V_PM_k = -1i*omega*conjugate(z_m - Pitching_Axis);
% Net Velocity
Vk = U + V_FV_k + V_PM_k;

%% Evaluating strength distribution of bound vortex sheet -----------------
% Applying kutta condition
[Str_Kutta,Circ_Kutta,Cp_Kutta,Fx_Kutta,Fy_Kutta] = Str_VS(z,Vk,U,0,"NO");
% Relaxed kutta condition
[K,R_Residual,R_Kutta] = K_fit(z,Vk,U,U + Vel_FV(Z2_centre,Circ,Radius,Z),Z,V);
[Str_Residual,Circ_Residual,Cp_Residual,Fx_Residual,Fy_Residual] = Str_VS(z,Vk,U,K,"NO");

%% Roots of strength distribution -----------------------------------------
% Arc Length
s = [0,cumsum(magnitude(z(2:end) - z(1:end-1)))];
% Surface Position - Functions
x_f = @(t) interp1(s,real(z),t);
y_f = @(t) interp1(s,imag(z),t);
% Evaluation with Kutta Condition
Str_Kutta_f = @(t) interp1(s,Str_Kutta,t);
s_Kutta = s(abs(diff(sign(Str_Kutta))) > 0);
z_Kutta = x_f(s_Kutta) + 1i*y_f(s_Kutta);
% Trailing Edge Limit with Kutta Condition
s_Kutta = [s(1),s_Kutta,s(end)];
z_Kutta = [z(1),z_Kutta,z(end)];
% Evaluation by Residual Minimisation
Str_Residual_f = @(t) interp1(s,Str_Residual,t);
s_Residual = s(abs(diff(sign(Str_Residual))) > 0);
z_Residual = x_f(s_Residual) + 1i*y_f(s_Residual);
% Trailing Edge Limit by Residual Minimisation
TE_LIMIT = min([abs(s_Residual),abs(s_Residual - s(end))]);

%% Normal Test ------------------------------------------------------------
% Note that non-zero values occur due to floating point errors
% Relative Vector Components
dx = real(z(2:end) - z(1:end-1));
dy = imag(z(2:end) - z(1:end-1));
% Inclination Angle
theta = atan2(dy,dx);
% Normal Vector Components (outwards)
Nx = -sin(theta);
Ny = cos(theta);
% Arc Length
s_m = 0.5*(s(2:end) + s(1:end-1));
% Normal components with Kutta Condition
Vk_Kutta = Vk(1:end-1) + Vel_VS(z,Str_Kutta,z_m(1:end-1));
uk_Kutta = real(Vk_Kutta)./sqrt(real(Vk_Kutta).*real(Vk_Kutta) + (-imag(Vk_Kutta)).*(-imag(Vk_Kutta)));
vk_Kutta = -imag(Vk_Kutta)./sqrt(real(Vk_Kutta).*real(Vk_Kutta) + (-imag(Vk_Kutta)).*(-imag(Vk_Kutta)));
N_Kutta = sum((uk_Kutta.*Nx + vk_Kutta.*Ny).*s_m);
% Normal components by Residual Minimisation
Vk_Residual = Vk(1:end-1) + Vel_VS(z,Str_Residual,z_m(1:end-1));
uk_Residual = real(Vk_Residual)./sqrt(real(Vk_Residual).*real(Vk_Residual) + (-imag(Vk_Residual)).*(-imag(Vk_Residual)));
vk_Residual = -imag(Vk_Residual)./sqrt(real(Vk_Residual).*real(Vk_Residual) + (-imag(Vk_Residual)).*(-imag(Vk_Residual)));
N_Residual = sum((uk_Residual.*Nx + vk_Residual.*Ny).*s_m);

%% Evaluating Flow Field --------------------------------------------------
% Free Vortices
V_FV = Vel_FV(Z2_centre,Circ,Radius,Z);
% Bound Vortex Sheet - Kutta Condition
V_BVS_Kutta = Vel_VS(z,Str_Kutta,Z);
% Bound Vortex Sheet - Residual Minimisation
V_BVS_Residual = Vel_VS(z,Str_Residual,Z);
% Flow Field - Kutta Condition
V_Kutta = U + V_FV + V_BVS_Kutta;
% Flow Field - Residual Minimisation
V_Residual = U + V_FV + V_BVS_Residual;

%% Setting Reference Frame ------------------------------------------------
% Experimental Data
V = V - 1i*omega*conjugate(Z - Pitching_Axis);
u = real(V);
v = -imag(V);
% Potential Flow Solution
V_Kutta = V_Kutta - 1i*omega*conjugate(Z - Pitching_Axis);
V_Residual = V_Residual - 1i*omega*conjugate(Z - Pitching_Axis);

%% Non-Dimensionalising Results -------------------------------------------
% Spatial
x = x/Chord_Length;
y = y/Chord_Length;
Z = Z/Chord_Length;
z = z/Chord_Length;
Pitching_Axis = Pitching_Axis/Chord_Length;
Z2_centre = Z2_centre/Chord_Length;
Radius = Radius/Chord_Length;
s = s/Chord_Length;
s_Kutta = s_Kutta/Chord_Length;
z_Kutta = z_Kutta/Chord_Length;
s_Residual = s_Residual/Chord_Length;
z_Residual = z_Residual/Chord_Length;
% Velocity
V = V/U;
u = u/U;
v = v/U;
Circ = Circ/(Chord_Length*U);
V_Kutta = V_Kutta/U;
V_Residual = V_Residual/U;
% Printed Results
Circ_Kutta = Circ_Kutta/(U*Chord_Length);
Circ_Residual = Circ_Residual/(U*Chord_Length);
K = K/U;
TE_LIMIT = TE_LIMIT/Chord_Length;

%% Residual Plot ----------------------------------------------------------
figure(1)
% Kutta Condition
subplot(1,2,1), shading interp, axis equal
surf(real(Z),imag(Z),abs((magnitude(V_Kutta) - magnitude(V))./magnitude(V)))
title('Relative Error - Kutta Condition')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
c = colorbar;
title(c,'\epsilon')
shading interp, view(2);
% Residual Minimisation
subplot(1,2,2), shading interp, axis equal
surf(real(Z),imag(Z),abs((magnitude(V_Residual) - magnitude(V))./magnitude(V)))
title('Relative Error - Residual Minimisation')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
c = colorbar;
title(c,'\epsilon')
shading interp, view(2);

%% Quiver/Streamline Plot -------------------------------------------------
figure(3), A = 10; B = 6; 
% Experimental Data
subplot(1,3,1), hold on, axis equal
plot(Z2_centre,'ro','MarkerFaceColor',[1,0,0])
quiver(x(1:B:end,1:B:end),y(1:B:end,1:B:end),u(1:B:end,1:B:end),v(1:B:end,1:B:end))
h = streamline(x,y,u,v,min(min(x))*ones(1,A),linspace(min(min(y)),max(max(y)),A));
fill(real(z),imag(z),'w')
set(h,'color','r')
title('Flow Field - Experimental Data')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
legend('Vortex Centre','location','southeast')
% Kutta Condition
subplot(1,3,2), hold on, axis equal
plot(Z2_centre,'ro','MarkerFaceColor',[1,0,0])
quiver(x(1:B:end,1:B:end),y(1:B:end,1:B:end), real(V_Kutta(1:B:end,1:B:end)),-imag(V_Kutta(1:B:end,1:B:end)))
h = streamline(x,y,real(V_Kutta),-imag(V_Kutta),min(min(x))*ones(1,A),linspace(min(min(y)),max(max(y)),A));
fill(real(z),imag(z),'w')
set(h,'color','r')
title('Flow Field - Kutta Condition')
legend('Vortex Centre','location','southeast')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
legend('Vortex Centre')
% Residual Minimisation
subplot(1,3,3), hold on, axis equal
plot(Z2_centre,'ro','MarkerFaceColor',[1,0,0])
quiver(x(1:B:end,1:B:end),y(1:B:end,1:B:end), real(V_Residual(1:B:end,1:B:end)),-imag(V_Residual(1:B:end,1:B:end)))
h = streamline(x,y,real(V_Residual),-imag(V_Residual),min(min(x))*ones(1,A),linspace(min(min(y)),max(max(y)),A));
fill(real(z),imag(z),'w')
set(h,'color','r')
title('Flow Field - Relaxed Kutta Condition')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
legend('Vortex Centre','location','southeast')

%% Velocity Magnitude Plot ------------------------------------------------
figure(4)
% Experimental Data
subplot(1,3,1), hold on, axis equal
surface(x,y,zeros(size(x)),sqrt(u.^2 + v.^2)), shading interp
fill3(real(z),imag(z),2*ones(size(z)),'w')
title('Velocity Magnitude - Experimental Data')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
c = colorbar;
title(c,'|v|/U')
caxis([0,max(max(sqrt(u.^2 + v.^2)))])
% Kutta Condition
subplot(1,3,2), hold on, axis equal
plot(z_Kutta,'ro','MarkerFaceColor','r')
surface(x,y,zeros(size(x)),sqrt(real(V_Kutta).^2 + (-imag(V_Kutta)).^2)), shading interp
fill(real(z),imag(z),'w')
plot(z_Kutta,'ro','MarkerFaceColor','r')
title('Velocity Magnitude - Kutta Condition')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
c = colorbar;
title(c,'|v|/U')
caxis([0,max(max(sqrt(u.^2 + v.^2)))])
legend('Separation/Attachment')
% Residual Minimisation
subplot(1,3,3), hold on, axis equal
plot(z_Residual,'ro','MarkerFaceColor','r')
surface(x,y,zeros(size(x)),sqrt(real(V_Residual).^2 + (-imag(V_Residual)).^2)), shading interp
fill(real(z),imag(z),'w')
plot(z_Residual,'ro','MarkerFaceColor','r')
title('Velocity Magnitude - Residual Minimisation')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
c = colorbar;
title(c,'|v|/U')
caxis([0,max(max(sqrt(u.^2 + v.^2)))])
legend('Separation/Attachment')

%% Gamma2 Criterion -------------------------------------------------------
figure(5), hold on, axis equal
surface(real(Z),imag(Z),zeros(size(Z)),Gamma2), shading interp
% Vortex Fit
for i = 1:length(Circ)
    text(real(Z2_centre(i)),imag(Z2_centre(i)),num2str(i),'Color','r','FontWeight','bold','HorizontalAlignment','Center')
    plot(real(Z2_centre(i)) + Radius(i)*cos(linspace(0,2*pi,100)),imag(Z2_centre(i)) + Radius(i)*sin(linspace(0,2*pi,100)),'r-','LineWidth',2)
end
fill(real(z),imag(z),'w')
title('Conditioned \gamma_2 - Criterion')
xlabel('x/C')
xlim([min(min(x)),max(max(x))])
xticks(linspace(min(min(x)),max(max(x)),5))
xticklabels(linspace(round(min(min(x)),1),round(max(max(x)),1),5))
ylabel('y/C')
ylim([min(min(y)),max(max(y))])
yticks(linspace(min(min(y)),max(max(y)),5))
yticklabels(linspace(round(min(min(y)),1),round(max(max(y)),1),5))
set(gca, 'YDir','reverse')
c = colorbar;
title(c,'\gamma_2')
caxis([0,1])

%% Bound Sheet Strength ---------------------------------------------------
figure(6)
subplot(2,1,1), hold on, grid on
plot(s,log(abs(Str_Kutta)),'b')
for i = 1:length(s_Kutta)
    plot([s_Kutta(i),s_Kutta(i)],1.2*[min(log(abs(Str_Kutta))),max(log(abs(Str_Kutta)))],'k:','LineWidth',2)
end
xlabel('s/C'), xlim([-0.2*s(end),1.2*s(end)])
ylabel('log(|\lambda|/U)'), ylim(1.2*[min(log(abs(Str_Kutta))),max(log(abs(Str_Kutta)))])
title('Bound Sheet Strength - Kutta Condition')
legend(['\Gamma^* = ',num2str(Circ_Kutta)],'Roots','location','southeast')
subplot(2,1,2), hold on, grid on
plot(s,log(abs(Str_Residual)),'r')
for i = 1:length(s_Residual)
    plot([s_Residual(i),s_Residual(i)],1.2*[min(log(abs(Str_Residual))),max(log(abs(Str_Residual)))],'k:','LineWidth',2)
end
xlabel('s/C'), xlim([-0.2*s(end),1.2*s(end)])
ylabel('log(|\lambda|/U)'), ylim(1.2*[min(log(abs(Str_Residual))),max(log(abs(Str_Residual)))])
title('Bound Sheet Strength - Relaxed Kutta Condition')
legend(['\Gamma^* = ',num2str(Circ_Residual)],'Roots','location','southeast')

%% Displaying Results -----------------------------------------------------
fprintf('\n')
fprintf('R Squared - Kutta Condition = %0.4f',R_Kutta)
fprintf('\n')
fprintf('Circulation - Kutta Condition = %0.4f [(m2/s)/(m2/s)]',Circ_Kutta)
fprintf('\n\n')
fprintf('R Squared - Residual Minimisation = %0.4f',R_Residual)
fprintf('\n')
fprintf('Circulation - Residual Minimisation = %0.4f [(m2/s)/(m2/s)]',Circ_Residual)
fprintf('\n\n')
fprintf('Trailing Edge Condition - K = %0.4f [(m/s)/(m/s)]',K)
fprintf('\n')
fprintf('Trailing Edge Limit - s = %0.4f [m/m]',TE_LIMIT)
fprintf('\n')
fprintf('Circulation Difference = %0.4f [(m2/s)/(m2/s)]',abs((Circ_Kutta - Circ_Residual)/(Circ_Kutta)))
fprintf('\n\n')
fprintf('Normal Velocity Integral - Kutta Condition = %0.4f',N_Kutta)
fprintf('\n')
fprintf('Normal Velocity Integral - Residual Minimisation = %0.4f',N_Residual)
fprintf('\n\n')
fprintf('Detection Window Size - Gamma1 = %0.4f [m/m]',L_1/Chord_Length)
fprintf('\n')
fprintf('Detection Window Size - Gamma2 = %0.4f [m/m]',L_2/Chord_Length)
fprintf('\n\n\n')
fprintf('Gamma1 Detection Ratios:')
disp(Ratio_1)
fprintf('Gamma2 Detection Ratios:')
disp(Ratio_2)
fprintf('\n')
fprintf('Circulation of Detected Vortices - [(m2/s)/(m2/s)] :')
disp(Circ)
fprintf('Relative Error in Circulation:')
disp(Circ_Error)
fprintf('\n')
fprintf('Radius of Detected Vortices - [m/m] :')
disp(Radius)
fprintf('Relative Error in Radius:')
disp(Radius_Error)
