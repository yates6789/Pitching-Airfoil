% Str_VS
% James Yates - 11/11/2019
% This script impliments the vortex panel method to evaluate the strength 
% distribution along a bound vortex sheet
% This is a modified version of the script written by C. Pozrikidis:
% Pozrikidis, C., 2017. Fluid Dynamics - Theory, Computation, and 
% Numerical Simulation. 3rd ed. s.l.:Springer.

%% Variables --------------------------------------------------------------
% z - Complex Coordinates of vortex sheet (must be closed loop) - [1xn]
% Vk - Complex Conjugate Velocity associated with z (without sheet) [1xn]
% U - Free Stream Velocity - [1]
% K - Net Strength of Trailing Edge (K = 0 for Kutta Condition) - [1]
% Setting PLOTS = 1 will produce plots

%% Outputs ----------------------------------------------------------------
% Str - Strength distribution of vortex sheet - [1xn]
% Circ - Circulation of vortex sheet - [1]
% Cp - Pressure Coefficient across airfoil surface - [1xn]
% Fx - Force per unit length in x-direction - [1]
% Fy - Force per unit length in y-direction - [1]

function [Str,Circ,Cp,Fx,Fy] = Str_VS(z,Vk,U,K,PLOTS)

%% Number of Vortex Panels ------------------------------------------------
PANELS = length(z) - 1;

%% Global Coordinates -----------------------------------------------------
X = real(z);
Y = imag(z);

%% Velocity Components ----------------------------------------------------
Ux = real(Vk); 
Uy = -imag(Vk);

%% Start and End Points of Panels -----------------------------------------
START = [X(1:end-1);Y(1:end-1)];
END = [X(2:end);Y(2:end)];

%% Surface Properties -----------------------------------------------------
% Relative Vector Components
dx = START(1,:) - END(1,:);
dy = START(2,:) - END(2,:);
dl = sqrt(dx.*dx + dy.*dy);
% Inclination Angle
theta = atan2(dy,dx);
% Tangent Vector Components
Sx = -cos(theta);
Sy = -sin(theta);
% Normal Vector Components (outwards)
Nx = sin(theta);
Ny = -cos(theta);

%% Collocation Points (panel mid-points) ----------------------------------
co = 0.5*(START+END);

%% Defining Global Influence Matrix ---------------------------------------
a = zeros(length(X),length(X));
b = zeros(length(X),length(X));
% Stepping through Collocation Points
for i = 1:PANELS
    % Stepping through Panels
    for j = 1:PANELS
        % Local Coordinates -----------------------------------------------
        % Collocation Point
        xt = co(1,i) - START(1,j);
        yt = co(2,i) - START(2,j);
        x = xt*Sx(j) + yt*Sy(j);
        y = -xt*Sy(j) + yt*Sx(j);
        % End Point
        x2t = END(1,j) - START(1,j);
        y2t = END(2,j) - START(2,j);
        x2 = x2t*Sx(j) + y2t*Sy(j);   
        % Subtended Values ------------------------------------------------
        % Radius Values
        r1 = sqrt(x*x+y*y);
        r2 = sqrt((x-x2)*(x-x2)+y*y);
        % Inclinaton Angles
        theta_1 = atan2(y,x);
        theta_2 = atan2(y,x-x2);
        % Local Influence Coefficients ------------------------------------
        if (i == j) % Correcting self-induced velocity                                                         
           ax1 = 0.5*(x/x2-1.0);
           ay1 = 1.0/(2*pi);
           ax2 = -0.5*x/x2;
           ay2 = -1.0/(2*pi);
        else
           dth = theta_2 - theta_1;
           ax1 = (1/(2*pi*x2))*( y*(log(r2/r1)) + (x-x2)*dth );
           ay1 = (1/(2*pi*x2))*((x-x2)*(log(r2/r1)) - y*dth + x2);
           ax2 = -(1/(2*pi*x2))*(y*(log(r2/r1)) + x*dth );
           ay2 = -(1/(2*pi*x2))*(x*(log(r2/r1)) - y*dth + x2);
        end     
       % Global Influence Coefficients ------------------------------------
       ux1 = ax1*Sx(j) - ay1*Sy(j);
       uy1 = ax1*Sy(j) + ay1*Sx(j);
       ux2 = ax2*Sx(j) - ay2*Sy(j);
       uy2 = ax2*Sy(j) + ay2*Sx(j);       
       % Velocity influence projected into normal vector ------------------
       if(j == 1) % Correcting self-induced velocity 
          a(i,1) = ux1*Nx(i) + uy1*Ny(i);
          Previous_Panel_a = ux2*Nx(i) + uy2*Ny(i);                      
       elseif(j == PANELS)
          a(i,PANELS) = ux1*Nx(i) + uy1*Ny(i) + Previous_Panel_a;
          a(i,PANELS+1) = ux2*Nx(i) + uy2*Ny(i);
       else
          a(i,j) = ux1*Nx(i) + uy1*Ny(i) + Previous_Panel_a;
          Previous_Panel_a = ux2*Nx(i) + uy2*Ny(i);                      
       end
       % Velocity influence projected into tangent vector -----------------
       if(j == 1) % Correcting self-induced velocity 
          b(i,1)= ux1*Sx(i) + uy1*Sy(i);
          Previous_Panel_b = ux2*Sx(i) + uy2*Sy(i);                      
       elseif(j == PANELS)
          b(i,PANELS) = ux1*Sx(i) + uy1*Sy(i) + Previous_Panel_b;
          b(i,PANELS+1) = ux2*Sx(i) + uy2*Sy(i);
       else
          b(i,j) = ux1*Sx(i) + uy1*Sy(i) + Previous_Panel_b;
          Previous_Panel_b = ux2*Sx(i) + uy2*Sy(i);                      
       end       
    end
end

%% Applying the Kutta Condition -------------------------------------------
% Both are positive direction vectors are in opposite direction 
a(PANELS+1,1) = 1;
a(PANELS+1,PANELS+1) = 1;

%% Normal Velocity Component ----------------------------------------------
% K is relaxation factor to allow tangential and normal flow
N_dot_U = [-Ux(1:end-1).*Nx - Uy(1:end-1).*Ny,K]';

%% Strength Calculation ---------------------------------------------------
Str = (a\N_dot_U)';

%% Pressure Coefficent Calculation ----------------------------------------
Tangent_Velocity = zeros(1,PANELS);
for i=1:PANELS
    Tangent_Velocity(i) = Ux(i)*Sx(i) + Uy(i)*Sy(i);
    for j = 1:PANELS+1
        Tangent_Velocity(i) = Tangent_Velocity(i) + b(i,j)*Str(j);
    end
end
Cp = 1.0 - (Tangent_Velocity.^2)/(U^2);
Cp = [Cp,Cp(1)];

%% Circulation Calculation ------------------------------------------------
Circ = 0.5*sum((Str(2:end) + Str(1:end-1)).*dl);

%% Force Coefficent Calculations - (F' = F/(0.5*rho*U^2) ------------------
Fx = sum(Cp(1:PANELS).*Nx.*dl);
Fy = sum(Cp(1:PANELS).*Ny.*dl);

%% Surface Vecotor Plots --------------------------------------------------
if PLOTS == 1
   figure(101), hold on, grid on, axis equal
   plot(X,Y,'k')
   quiver(co(1,:),co(2,:),Nx,Ny)
   quiver(co(1,:),co(2,:),Sx,Sy)
   title('Surface Vectors')
   xlabel('x - [m]'), xlim([min(X) - 0.5*(max(X) - min(X)),max(X) + 0.5*(max(X) - min(X))])
   ylabel('y - [m]'), ylim([min(Y) - 0.5*(max(X) - min(X)),max(Y) + 0.5*(max(X) - min(X))])
   legend('Surface','Normal Vector','Tangent Vector')
end

%% Strength Plot ----------------------------------------------------------
if PLOTS == 1
   figure(102), hold on
   plot([0,cumsum(sqrt(dx(1+m:end-m).^2 + dy(1:end).^2))],Str(1:end))
   title('Panel Strength')
   xlabel('s - [m]'), xlim([0,max(max(cumsum(sqrt(dx.*dx + dy.*dy))))])
   ylabel('\gamma - [m/s]'), ylim(1.1*[min(Str(1:end)),max(Str(1:end))])
   grid on
end

%% Cp Plot ----------------------------------------------------------------
if PLOTS == 1
   figure(103), hold on
   plot(X(1+m:end-m),Cp(1+m:end-m))
   title('Pressure Coefficient')
   xlabel('x - [m]'), xlim([min(X) - 0.1*(max(X)),1.1*max(X)])
   ylabel('Cp - []'), ylim([min(Cp) - 0.1*(max(Cp)),1.1*max(Cp)])
   grid on
   set(gca,'ydir','reverse')
end
end
