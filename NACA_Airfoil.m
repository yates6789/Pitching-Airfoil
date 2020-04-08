% NACA_Airfoil
% James Yates - 08/02/2020
% This script generates the complex coorindates of a NACA airfoils
% This is a modified version of the script written by C. Pozrikidis:
% Pozrikidis, C., 2017. Fluid Dynamics - Theory, Computation, and 
% Numerical Simulation. 3rd ed. s.l.:Springer.

%% Variables --------------------------------------------------------------
% Maximum_Camber - Maximum camber of the airfoil - [1]
% Maximum_Camber_Position - Position of Maximum Camber - [1]
% thickenss - Thickness of the airfoil - [1]
% Chord_Length - Chord Length of the Airfoil - [1]
% Pitching_Axis - Position of pitching axis in complex coordinates - [1]
% theta - Angle of attack - [1] - [rad]
% PANELS - Number of Panels for discretisation - [1]
% Setting PLOTS = "YES" will produce plots

%% Outputs ----------------------------------------------------------------
% za - complex coordinates of the airfoil surface - [1x(PANELS+1)]

function za = NACA_Airfoil(maximum_camber,maximum_camber_position,thickness,Chord_Length,Pitching_Axis,theta,PANELS,PLOTS)

%% Initialising coordinate vectors ----------------------------------------
% Airfoil Perimeter
xa = zeros(1,PANELS);
ya = zeros(1,PANELS);
% Camber Line 
yc = zeros(1,PANELS);

%% Stepping through discretisation
for i=1:PANELS
    psi = (i-1.0)*2*pi/PANELS;
    x = 0.5*(1.0+cos(psi));
    % Camber Line
    if(x<maximum_camber_position)
       yc(i) = maximum_camber*(2.0*maximum_camber_position*x-x*x)/maximum_camber_position^2;
    else
       yc(i) = maximum_camber/(1.0-maximum_camber_position)^2*((1.0-2.0*maximum_camber_position)+2.0*maximum_camber_position*x-x*x);
    end
    % Thickness Approximation
    yt = 5.0*thickness *(0.2969*sqrt(x)-0.1260*x-0.3516*x*x+0.2843*x*x*x-0.1036*x*x*x*x);
    % Profile Perimeter
    if(i<(PANELS/2+1))
       ya(i) = (yc(i)-yt)*Chord_Length;
    else
       ya(i) = (yc(i)+yt)*Chord_Length;
    end
    xa(i) = x*Chord_Length;
end

%% Closing Loop with Trailing Edge ----------------------------------------
xa(PANELS+1) = Chord_Length;
ya(PANELS+1) = 0.0;
yc(PANELS+1) = 0.0;

%% Centering around Pitching Axis -----------------------------------------
xa = xa - 0.25*Chord_Length;

%% Changing Pitch ---------------------------------------------------------
% Chord
[tc,rc] = cart2pol(xa,Chord_Length*yc);
[~,yc] = pol2cart(tc - theta,rc);
% Airfoil
[ta,ra] = cart2pol(xa,ya);
[xa,ya] = pol2cart(ta - theta,ra);

%% Centring around point of interest --------------------------------------
xa = xa + real(Pitching_Axis);
ya = ya + imag(Pitching_Axis);
yc = yc + imag(Pitching_Axis);

%% Complex Coordinates ----------------------------------------------------
za = xa + 1i*ya;

%% Plotting Results -------------------------------------------------------
if PLOTS == "YES"
   figure(201), hold on, axis equal
   plot(xa,ya,'k-');
   plot(xa,yc,'k--');
   xlabel('x')
   ylabel('y')
   grid on
end
                   
end
