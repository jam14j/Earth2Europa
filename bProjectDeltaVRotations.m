clear all 
close all
clc; 

vEaJ = [8.8867;0.3226;0.2017]; %SCI velocity arriving at Jupiter
vJEu = [-14.8552;-31.8424;0]; %synodic velocity leaving jupiter
muSun = 1.3271244004193938e11;
muJ=126.687e6; %km^3/s^2
rJ=71492; %km

arrivalDay = calToJD(2030,7,4,12,0,0); 
[rJArr,vJArr] = planetToPos(5,2030,7,4,12,0,0,muSun);

%Convert heliocentric vEarthJ to jupiter centric JCI frame
vEaJ_relJ = vEaJ-vJArr; 
jTilt = 3.13; %degrees
% %rotTilt = [1,0,0;
%     0,cos(jTilt),sin(jTilt);
%     0,-sin(jTilt),cos(jTilt)];
vEaJ_jci = vEaJ_relJ;%rotTilt*vEaJ_relJ; 

%Convert synodic vJEu to jci frame
ws = deg2rad(101.38)/(24*3600); %rad per second
ti = 0; 
tstep = 100; 
tf = 10000; 

day1 = calToJD(2013,10,10,4,30,0); 

diffInDays = arrivalDay-day1; 
arrivalAngle = ws*diffInDays; 
vJEu_jci = rotz(rad2deg(arrivalAngle))*vJEu + cross([0;0;ws],vJEu); 

% %Calculate the required delta v
% rpInit = [-100000,200000];%149236.008237;
% options = optimset('TolX',1e-12); 
% dTotal = acos(dot(vEaJ_jci,vJEu_jci)/(norm(vEaJ_jci)*norm(vJEu_jci))); 
% [rpFin] = fzero(@(rp)f(rp,dTotal,muJ,vEaJ_jci,vJEu_jci),rpInit,options); 
% disp(rpFin)
% 
% vbneg = sqrt(norm(vEaJ_jci)^2+2*muJ/rpFin);
% vbpos = sqrt(norm(vJEu_jci)^2+2*muJ/rpFin);
% dV = abs(vbpos-vbneg);
% disp(dV)


rp = 149236.008237;
%Convert heliocentric vEarthJ to jupiter centric JCI frame
r = [rp,0,0]; 
state = [r,transpose(vJEu_jci)]; 
tspan = -3600:10:36000;
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, x] = ode113(@(t1,x1) FODE(t1,x1,muJ), tspan, state);
[val,ind] = min(vecnorm(transpose(x(:,1:3)))');
r2 = x(ind,1:3);%[val,0,0];
state2 = [r2,transpose(vEaJ_jci)]; 
tspan2 = -3600:10:36000;
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t2, x2] = ode113(@(t1,x1) FODE(t1,x1,muJ), tspan2, state2);
figure
circle(0,0,rJ)
hold on
plot3(x(:,1),x(:,2),x(:,3));
hold on
plot3(x2(:,1),x2(:,2),x2(:,3));

function jd = calToJD(y,m,d,h,min,s)
%Converst from date/time to JD
jd = datenum(y,m,d,h,min,s)+1721058.5;
end

function [r,v] = planetToPos(planet,y,m,d,h,min,s,mu)
%Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
%Code takes in planet number and time/date and converts to r and v

%Get elements from given code
[elements,rates]= AA279j2000_planetary_elements(planet);

%Get elements out from data input
a = elements(1);
e = elements(2);
i = elements(3);
Om = elements(4);
w_bar = elements(5);
L = elements(6);
adot = rates(1);
edot = rates(2);
idot = rates(3);
Omdot = rates(4);
w_bardot = rates(5);
Ldot = rates(6);

jd = calToJD(y,m,d,h,min,s); %Calculate JD from input time/date
%Convert J2000 planetary orbits/rates to classic orbit elements
a_t = a+adot*(jd-2451545)/36525; 
e_t = e+edot*(jd-2451545)/36525; 
i_t = deg2rad(i+idot*(jd-2451545)/36525); 
Om_t = deg2rad(Om+Omdot*(jd-2451545)/36525); 
w_bar_t = w_bar+w_bardot*(jd-2451545)/36525; 
L_t = L + Ldot*(jd-2451545)/36525; 

w_t = deg2rad(w_bar_t) - Om_t; 
M_t = deg2rad(L_t-w_bar_t); 

E = meanToEcc(M_t,e,1e-9); %convert mean to eccentric anomaly

 %convert to km and return values
[r,v] = oeToSCI(a_t*149597870.7,e_t,i_t,Om_t,w_t,E,mu);
end

function [r_out,v_out] = oeToSCI(a,e,i,Om,w,E,mu)
%Converts classic orbital elements to r and v in sun centered inertial
%frame
rpqw=[a*(cos(E)-e);
    a*sqrt(1-e^2)*sin(E);
    0];
vpqw=[a*sqrt(mu/(a^3))/(1-e*cos(E))*(-sin(E)); 
   a*sqrt(mu/(a^3))/(1-e*cos(E))*(sqrt(1-e^2)*cos(E));
   0];
Rpqw_ijk=rotz(rad2deg(Om))*rotx(rad2deg(i))*rotz(rad2deg(w)); 
r_out=(Rpqw_ijk*rpqw); 
v_out=(Rpqw_ijk*vpqw);
end

function [year,mon,day,hour,min,sec,cal] = JDtoCal(jd)
%Converts from JD to date/time
calStr = datestr(jd-1721058.5);
year = str2double(calStr(8:11)); 
day = str2double(calStr(1:2)); 
hour = str2double(calStr(13:14));
min = str2double(calStr(16:17)); 
sec = str2double(calStr(19:20)); 
mon = month(datetime(calStr(1:11))); 
cal = [year,mon,day,hour,min,sec];
end

function dx = FODE(t, x, mu)
    dx = zeros(6,1);
    dx(1:3) = x(4:6);
    dx(4:6) = -mu*x(1:3)/norm(x(1:3))^3;
end

%Function to go to zero to solve for rp
function [y] = f(rp,dTot,mu,v1,v2)
y = asin(1/(1+rp*norm(v1)^2/mu))+asin(1/(1+rp*norm(v2)^2/mu))-dTot;
if ~isreal(y)
    disp('hi')
end
end