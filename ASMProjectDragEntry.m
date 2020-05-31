% AA279B - Class Project
% From Earth to the Water on Europa
% Part 3: Europa Landing
% Code by:
% Christine Hamilton and
% Juan Martinez Castellanos
clc; clear all; close all;
% Constants
Reu = 3122000/2; % [m]
Meu = 4.799844*10^22; % [kg]
G = 6.67430*10^-11; % [m3/kg1s2]
mueu = G*Meu;
Msat = 600; % [kg]

% Simulation
r0 = 65000+Reu; % [m]
v0 = sqrt(G*Meu/r0);
x0 = [r0;0;0;0;v0;0];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t,xout] = ode113(@(t,x) FODEwDrag(t,x,G*Meu,Reu,Msat), 0:100:100*24*3600, x0, opts);
rout = vecnorm(xout(:,1:3)')';
% Plot Orbit
figure
hold on
plot(xout(:,1)/1000,xout(:,2)/1000);
% Plot circle
th = 0:pi/50:2*pi;
xcirc = Reu/1000*cos(th);
ycirc = Reu/1000*sin(th);
plot(xcirc, ycirc);
hold off
axis equal
title("Reentry Using Drag");
xlabel("x (km)")
ylabel("y (km)")
legend("Trajectory", "Europa"); 


% FUNCTIONS
function dx = FODEwDrag(t, x, mu, Reu, Msat)
    r = x(1:3); v = x(4:6);
    g = -mu/norm(r)^2;
    drag = bEuropaAtmosphericModel(norm(r)-Reu,norm(v),g);
    if ~isreal(drag)
        fprintf("COMPLEX DRAG\n");
    end
%     drag = 0;
    dx = zeros(6,1);
    dx(1:3) = v;
    dx(4:6) = -mu*r/norm(r)^3 - (drag/Msat)*(v/norm(v));
end



