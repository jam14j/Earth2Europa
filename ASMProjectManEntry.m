% AA279B - Class Project
% From Earth to the Water on Europa
% Part 3: Europa Landing with Maneuver
% Code by:
% Christine Hamilton and
% Juan Martinez Castellanos
clc; clear all; close all;
% Constants
Reu = 3122000/2; % [m]
Meu = 4.799844*10^22; % [kg]
G = 6.67430*10^-11; % [m3/kg1s2]
mueu = G*Meu;

% Initial conditions
r0 = 300000+Reu; % [m]
r1 = [r0;0;0];
v0 = sqrt(G*Meu/r0);
v0 = [0;v0;0];
% Target
r2 = [0;0;Reu];

% Lambert to surface
% find minimum time
% c = norm(rf_ECI-r0_ECI);
% s = (norm(rf_ECI)+norm(r0_ECI))/2;
% t_min = sqrt(2/)
collision_flag = 0;
dV_opt = 99999999999;
for t = 10:10:10000
    [v1, v2] = AA279lambert_curtis(mueu,r1,r2,'pro',0,t);
    % Check for collision
    collision_flag = 0;
%     if ((dot(r1,v1)<0) && (dot(r2,v2)>0))
%         energy = norm(v1)^2/2-mueu/norm(r1);
%         a = -mueu/(2*energy);
%         h = norm(cross(r1,v1));
%         p = h^2/mueu;
%         e = sqrt((a-p)/a);
%         rp = a*(1-e);
%         if rp <= Reu
%             collision_flag = 1;
%         end
%     end
    % TEST TEST
    x0 = [r1;v1];
    [ttest, xtest] = ode113(@(t1,x1) FODE(t1,x1,mueu), 0:10:t, x0);
    rp = min(vecnorm(xtest(:,1:3)')');
    if rp <= Reu
        collision_flag = 1;
    end
    
    
    dV_tot = norm(v1-v0) + norm(v2);
    if (collision_flag==0) && (dV_tot<dV_opt)
        v1_opt = v1;
        v2_opt = v2;
        dV_opt = dV_tot;
        t_opt = t;
    end
end
fprintf("Optimal dV = %f km/sec\n", dV_opt/1000);
fprintf("Optimal time = %f sec\n", t_opt);

% Simulation
x0 = [r1;v0];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, x] = ode113(@(t1,x1) FODE(t1,x1,mueu), 0:10:10000, x0);
x0 = [r1;v1_opt];
[tman, xman] = ode113(@(t1,x1) FODE(t1,x1,mueu), 0:10:t_opt, x0);

%% Plots
[xeu , yeu, zeu] = ellipsoid(0, 0, 0, Reu, Reu, Reu, 20);
figure
hold on
line(x(:,1)/1000, x(:,2)/1000, x(:,3)/1000, 'Color','r');
line(xman(:,1)/1000, xman(:,2)/1000, xman(:,3)/1000,'Color','g');
surface(xeu/1000, yeu/1000, zeu/1000, 'EdgeColor', 'black');
hold off
xlabel("x (km)")
ylabel("y (km)")
zlabel("z (km)")
title("Trajectory of the Probe");
legend("Parking orbit","Landing trajectory")
axis equal;
    
% FUNCTIONS
function dx = FODE(t, x, mu)
    dx = zeros(6,1);
    dx(1:3) = x(4:6);
    dx(4:6) = -mu*x(1:3)/norm(x(1:3))^3;
end

function plotAroundEarth(r, tl)
    Re = 6.3781*10^6; % [m]
    [xE , yE, zE] = ellipsoid(0, 0, 0, Re, Re, Re, 20);
    figure
    hold on
    surface(xE, yE, zE, 'EdgeColor', 'black');
    line(r(:,1), r(:,2), r(:,3));
    xlabel("x (km)")
    ylabel("y (km)")
    zlabel("z (km)")
    title(tl)
    axis equal;
end

