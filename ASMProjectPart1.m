% AA279B - Class Project
% From Earth to the Water on Europa
% Part 1: Optimal Earth-Jupiter transfer orbit
% Code by:
% Christine Hamilton and
% Juan Martinez Castellanos
clc; clear all; close all;
% Constants
AU = 149597870.7; % [km]
muS =  1.3271244004193938*10^11; % [km3/s2]

% Launch and Arrival Times
% yyyy,mm,dd,hh,mm,ss
JD_initial = datenum(2025,1,1) + 1721058.5;
LTs = JD_initial:10:JD_initial+365.25;
ATs = JD_initial+(365.25*2):10:JD_initial+(365.25*8);

% Calculate delta-Vs
dvs = zeros(length(LTs), length(ATs));
vfs = zeros(length(LTs), length(ATs),3);
v0s = zeros(length(LTs), length(ATs),3);
for Lidx = 1:length(LTs)
    [rE, vE] = findPlanet(3,LTs(Lidx));
    for Aidx = 1:length(ATs)
        [rJ, vJ] = findPlanet(5,ATs(Aidx));
        [v0, vf] = AA279lambert_curtis(muS,rE,rJ,'pro', 0,...
                                      (ATs(Aidx)-LTs(Lidx))*24*3600);
        dvs(Lidx,Aidx) = norm(v0-vE);% + norm(vf-vJ);
        v0s(Lidx,Aidx,:) = v0;
        vfs(Lidx,Aidx,:) = vf;
    end
end

% Recover dates
LTs = LTs - 1721058.5;
ATs = ATs - 1721058.5;

% Minimize
[minrow, minj] = min(dvs);
[mindv, mink] = min(minrow);
optAT = ATs(mink);
optLT = LTs(minj(mink));
optvf = zeros(3,1);
optvf(:) = vfs(minj(mink),mink,:);
optv0 = zeros(3,1);
optv0(:) = v0s(minj(mink),mink,:);

fprintf("Minimum deltaV = %f km/s\n", mindv);
fprintf("Opt launch date = %s\n", datestr(optLT, 'mm/DD/YYYY'));
fprintf("Opt arrival date = %s\n", datestr(optAT, 'mm/DD/YYYY'));
fprintf("Time of flight = %f years\n", (optAT-optLT)/365.25);
fprintf("Velocity on arrival = %f km/s\n", norm(optvf));
fprintf("Departure velocity = %f km/s\n", norm(optv0));


% Simulate orbit
[r0, v0] = findPlanet(3,optLT+1721058.5);
x0 = [r0; optv0];
tspan = 0:1000000:(optAT-optLT)*24*3600;
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, x] = ode113(@(t1,x1) FODE(t1,x1,muS), tspan, x0);

% Planet orbits
rE_list = zeros(length(t),3);
rJ_list = zeros(length(t),3);
for j = 1:length(t)
    [rE_list(j,:), v0] = findPlanet(3,optLT+1721058.5+t(j)/(24*3600));
    [rJ_list(j,:), v0] = findPlanet(5,optLT+1721058.5+t(j)/(24*3600));
end

% Plot orbit
[rE,vE] = findPlanet(3, optLT+1721058.5);
[rJ,vJ] = findPlanet(5, optAT+1721058.5);
figure
hold on
line(x(:,1), x(:,2));
plot(rE_list(:,1),rE_list(:,2), 'b');
plot(rJ_list(:,1),rJ_list(:,2), 'r');
% plot(rE(1),rE(2),'bo')
% plot(rJ(1),rJ(2),'ro')
plot(0,0,'y*')
hold off
set(gca,'Color','k')
xlabel("x (km)")
ylabel("y (km)")
lgd = legend("Path","Earth","Jupiter","Sun");
lgd.TextColor = 'white';
lgd.Location = 'southeast';
title("Minimum \DeltaV Transfer Orbit")
axis equal
ylim([-8e8,3e8])

% Plot contour
figure
hold on
contourf(ATs, LTs, dvs, 'ShowText', 'on')
plot(ATs(mink), LTs(minj(mink)), 'r.', 'MarkerSize', 20)
hold off
datetick('x', 'yyyy', 'keepticks')
datetick('y', 1, 'keepticks')
xlabel("Arrival Date")
ylabel("Launch Date")
title("Contour of \DeltaV")

figure
hold on
surf(ATs, LTs, dvs, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
plot3(ATs(mink), LTs(minj(mink)), mindv, 'r.', 'MarkerSize', 50)
hold off
datetick('x', 'yyyy', 'keepticks')
datetick('y', 1, 'keepticks')
title("Surface of \DeltaV")
xlabel("Arrival Date")
ylabel("Launch Date")
zlabel("\DeltaV (km/s)")

function [rSCI,vSCI] = findPlanet(planetID, JD)
    muS =  1.3271244004193938*10^11; % [km3/s2]
    J2000_JD = 2451545.0;
    AU =  149597870.7; %[km]
    t = JD-J2000_JD; % [days]
    [oe,doe] = AA279j2000_planetary_elements(planetID);
    a = (oe(1)+doe(1)*t/36525)*AU;   % [km]
    e = oe(2)+doe(2)*t/36525;
    i = oe(3)+doe(3)*t/36525;        % [deg]
    RAAN = oe(4)+doe(4)*t/36525;     % [deg]
    om_tilde = oe(5)+doe(5)*t/36525; % [deg]
    om = om_tilde - RAAN;            % [deg]
    M = deg2rad(oe(6)+doe(6)*t/36525 - om_tilde); % [rad]
    
    n = sqrt(muS/a^3);
    nu = KeplerSolver2Nu(M, e);
    [rSCI,vSCI] = OEtoSCI(n,e,i,RAAN,om,nu);
end

function nu = KeplerSolver2Nu(M, e)
    eps = 10^-10;
    E = pi;
    d = -(E-e*sin(E)-M)/(1-e*cos(E));
    while(abs(d)>eps)
        d = -(E-e*sin(E)-M)/(1-e*cos(E));
        E = E + d;
    end
    nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end

function [rSCI,vSCI] = OEtoSCI(n,e,i,RAAN,om,nu)
    muS =  1.3271244004193938*10^11; % [km3/s2]
    a = (muS/(n^2))^(1/3);
    P = a*(1-e^2);
    r = P/(1+e*cos(nu));
    rPQW = r*[cos(nu); sin(nu); 0];
    vPQW = sqrt(muS/P)*[-sin(nu); e+cos(nu); 0];
    rSCI = rotz(RAAN)*rotx(i)*rotz(om)*rPQW;
    vSCI = rotz(RAAN)*rotx(i)*rotz(om)*vPQW;
end

function dx = FODE(t, x, mu)
    dx = zeros(6,1);
    dx(1:3) = x(4:6);
    dx(4:6) = -mu*x(1:3)/norm(x(1:3))^3;
end