close all; 
clear all; 
clc; 

G=6.673e-20; %km^3/(kg*s^2)
massE=4.799844e22; %kg
muJ=126.687e6; %km^3/s^2
muE=G*massE; %km^3/s^2
rJ=71492; %km
rE=1561; %km
rJE=670900;%149236.008237 km/sec; %km


t=[3*24*3600]; %seconds
%%
%Part A
r1=[-muE/(muE+muJ)*(norm(rJE));0;0];
r2=[muJ/(muE+muJ)*(norm(rJE));0;0];
L1y=0;
x0=norm(r1);
L1x=fzero(@f,x0);
fprintf('\tPosition of L1: %fi km, %gj km\n', L1x, L1y)
L1x = 6.6988e5; 

%%
%Set initial LJO parameters
v_init=35*sqrt(149236.008237/1.4433e+06);%[-10:1:0];%0;%[0:5:50];%35%55%10.92367104; %km/s
dv=0.001; 
dphi=0.001; 
phi_init=0;%[0:pi/2dr:2*pi];%pi/2 %3*pi/2; %deg2rad(47.70061087); %rad
rInit= 1.4433e+06;%149236.008237;%+rJ;


R1=muE*rJE/(muJ+muE);
R2=muJ*rJE/(muJ+muE);
ws=sqrt((muJ+muE)/(R1+R2)^3);
rI1=[-R1;0;0];
rI2=[R2;0;0];
tol=5; %km
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
tspan=[0:10:t];

J=zeros(2,2); 
dr=100;
while norm(dr)>tol
% loop = 1;
% for j = v_init
%     for k = phi_init
%         High velocity
%         v_high=j+dv; 
%         rvH=[-rInit*cos((k))-R1;
%             -rInit*sin(k);
%             0];
%         vvH=[v_high*sin(k);
%             -v_high*cos(k);
%             0];
%         [tvH,yvH]=ode113(@dVectJup, tspan,[rvH,vvH],options);
%         len=length(tvH);
%         rv_high=[yvH(len,1),yvH(len,2)];
% 
%         low velocity
%         v_low=j-dv; 
%         rvL=[-rInit*cos((k))-R1;
%             -rInit*sin(k);
%             0];
%         vvL=[v_low*sin(k);
%             -v_low*cos(k);
%             0];
%         [tvL,yvL]=ode113(@dVectJup, [0:10:t],[rvL,vvL],options);
%         len2=length(tvL);
%         rv_low=[yvL(len2,1),yvL(len2,2)];
%         error=transpose(((rv_high-rv_low)/(2*dv)));
%         J(1,1)=error(1); 
%         J(2,1)=error(2);
% 
%         High phi
%         phi_high=k+dphi;
%         rphiH=[-rInit*cos((phi_high))-R1;
%             -rInit*sin(phi_high);
%             0];
%         vphiH=[j*sin(phi_high);
%             -j*cos(phi_high);
%             0];
%         [tphiH,yphiH]=ode113(@dVectJup, tspan,[rphiH,vphiH],options);
%         len=length(tphiH);
%         rphi_high=[yphiH(len,1),yphiH(len,2)];
% 
%         low phi
%         phi_low=k-dphi; 
%         rphiL=[-rInit*cos((phi_low))-R1;
%             -rInit*sin(phi_low);
%             0];
%         vphiL=[j*sin(phi_low);
%             -j*cos(phi_low);
%             0];
%         [tphiL,yphiL]=ode113(@dVectJup, [0:10:t],[rphiL,vphiL],options);
%         len2=length(tphiL);
%         rphi_low=[yphiL(len2,1),yphiL(len2,2)];
%         error=transpose(((rphi_high-rphi_low)/(2*dphi)));
%         J(1,2)=error(1); 
%         J(2,2)=error(2);
% 
%         r=[-rInit*cos((k))-R1;
%             -rInit*sin((k));
%             0];
%         v=[j*sin((k));
%             -j*cos((k));
%             0];
%         [tav,yav]=ode113(@dVectJup, [0:10:t],[r,v],options);
%         len3=length(tav);
%         rb_N=[yav(len3,1);yav(len3,2)];
%         dr=rb_N-[L1x;L1y]; 
%         mat=[j;k]-inv(J)*(dr);
%         v_init=mat(1); 
%         phi_init=mat(2);
%         fprintf('v_init = %f \n phi_init = %g\n', j,k)
%         dr_mat(loop,:) = [j,k,norm(dr)];
%         disp(norm(dr))
%         loop = loop+1;
%     end
% end

%disp(dr_mat)
    %High velocity
    disp(1)
    v_high=v_init+dv; 
    rvH=[-rInit*cos((phi_init))-R1;
        -rInit*sin(phi_init);
        0];
    vvH=[v_high*sin(phi_init);
        -v_high*cos(phi_init);
        0];
    [tvH,yvH]=ode113(@dVectJup, tspan,[rvH,vvH],options);
    len=length(tvH);
    rv_high=[yvH(len,1),yvH(len,2)];
    
    %low velocity
    disp(2)
    v_low=v_init-dv; 
    rvL=[-rInit*cos((phi_init))-R1;
        -rInit*sin(phi_init);
        0];
    vvL=[v_low*sin(phi_init);
        -v_low*cos(phi_init);
        0];
    [tvL,yvL]=ode113(@dVectJup, [0:10:t],[rvL,vvL],options);
    len2=length(tvL);
    rv_low=[yvL(len2,1),yvL(len2,2)];
    error=transpose(((rv_high-rv_low)/(2*dv)));
    J(1,1)=error(1); 
    J(2,1)=error(2);
    
    %High phi
    disp(3)
    phi_high=phi_init+dphi;
    rphiH=[-rInit*cos((phi_high))-R1;
        -rInit*sin(phi_high);
        0];
    vphiH=[v_init*sin(phi_high);
        -v_init*cos(phi_high);
        0];
    [tphiH,yphiH]=ode113(@dVectJup, tspan,[rphiH,vphiH],options);
    len=length(tphiH);
    rphi_high=[yphiH(len,1),yphiH(len,2)];
    
    %low phi
    disp(4)
    phi_low=phi_init-dphi; 
    rphiL=[-rInit*cos((phi_low))-R1;
        -rInit*sin(phi_low);
        0];
    vphiL=[v_init*sin(phi_low);
        -v_init*cos(phi_low);
        0];
    [tphiL,yphiL]=ode113(@dVectJup, [0:10:t],[rphiL,vphiL],options);
    len2=length(tphiL);
    rphi_low=[yphiL(len2,1),yphiL(len2,2)];
    error=transpose(((rphi_high-rphi_low)/(2*dphi)));
    J(1,2)=error(1); 
    J(2,2)=error(2);
    
    disp(5)
    r=[-rInit*cos((phi_init))-R1;
        -rInit*sin((phi_init));
        0];
    v=[v_init*sin((phi_init));
        -v_init*cos((phi_init));
        0];
    [tav,yav]=ode113(@dVectJup, [0:10:t],[r,v],options);
    len3=length(tav);
    rb_N=[yav(len3,1);yav(len3,2)];
    dr=rb_N-[L1x;L1y]; 
    mat=[v_init;phi_init]-inv(J)*(dr);
    v_init=mat(1); 
    phi_init=mat(2);
    disp(6)
    disp(norm(dr))
end
fprintf('\tInitial velocity %f km/s\n', v_init)
fprintf('\tInitial angle %f rad\n', phi_init)

r=[-rInit*cos((phi_init))-R1;
    -rInit*sin((phi_init));
    0];
v=[v_init*sin((phi_init));
    -v_init*cos((phi_init));
    0];
[t,y]=ode113(@dVectJup, [0:10:t],[r,v],options);
%%
figure
circle(-R1,0,rJ)
hold on
circle(R2,0,rE)
hold on
%circle(R2,0,L1x-rJE)
%plot(R2,0,'*')
hold on
%circle(-R1,0,149236)
%plot(y(:,1),y(:,2));
legend('Jupiter', 'Europa','Trajectory','Location','southeast')
title('Jupiter/Europa Synodic Frame')
xlabel('km')
ylabel('km')

vx=(y(length(t),4));
vy=(y(length(t),5)); 
vz=(y(length(t),6)); 
dV=vecnorm(y(length(t),4:6)) - 1.312028402;%(norm([-vx,1.312028402-vy,-vz]));
fprintf('\tThe scalar value of the retro-burn at L1 is %f km/s\n',dV); 
%%
dvTot=(v_init-sqrt(muJ/(muJ+149236-rJ)))+dV; 
fprintf('\tThe total propellant is %f km/s\n', dvTot);


rI3_tot=[];
vI3_tot=[];
for j=1:length(t)
    c_si=[cos(ws*t(j)),sin(ws*t(j)),0;-sin(ws*t(j)),cos(ws*t(j)),0;0,0,1];
    c_is=inv(rotz(-rad2deg(ws*t(j))));
    r13_s=[y(j,1); y(j,2); y(j,3)];
    rI3_i=c_is*r13_s;
    rI3_tot=[rI3_tot; transpose(rI3_i)]; 
    vI3_s=[y(j,4); y(j,5); y(j,6)]; 
    vI3_i=vI3_s+cross([0;0;ws], rI3_i);
    vI3_tot=[vI3_tot;transpose(vI3_i)];
end
%%
r_init_sun = [18505152.613661; -44398019.742906; 3070992.384946];% [km]
v_init_sun = [-2.952874; 6.984550; -0.481079]; 
tspan = [0:10:100*24*3600];

[tsun,ysun]=ode113(@dVectJCI, tspan,[r_init_sun,v_init_sun],options);

figure
plot(rI3_tot(1:1500,1),rI3_tot(1:1500,2));
hold on
circle(-R1,0,rJ);
hold on
plot(ysun(:,1),ysun(:,2))

[rMin,ind] = min(vecnorm(ysun(:,1:3)')');

%%
figure
plot(rI3_tot(:,1), rI3_tot(:,2)); 
hold on
axis equal
title('Jupiter Europa Trajectory in Inertial Frame')
xlabel('km')
ylabel('km')
circle(-R1,0,rJ)
legend('Trajectory','Jupiter','Location','southeast')

%%
a = rInit; 
e = 0; 
i = 0; 
Om = 0; 
w = 0; 
E = 0; 
[r,v] = oeToSCI(a,e,i,Om,w,E,muJ)

function y=f(x)
G=6.673e-20; %km^3/(kg*s^2)
massE=4.799844e22; %kg
muJ=126.687e6; %km^3/s^2
muE=G*massE; %km^3/s^2
rJE=670900; %km

r1=muE/(muE+muJ)*(rJE);
r2=muJ/(muE+muJ)*(rJE);
ws=sqrt((muJ+muE)/(r1+r2)^3);
y=ws^2*x+muE/(x-r2)^2-muJ/(x+r1)^2; 
end

function[statedot] = dVectJup(t,state)
statedot=zeros(6,1); 
G=6.673e-20; %km^3/(kg*s^2)
massE=4.799844e22; %kg
muJ=126.687e6; %km^3/s^2
muE=G*massE; %km^3/s^2
rJE=670900; %km

R1=muE*rJE/(muJ+muE);
R2=muJ*rJE/(muJ+muE);
ws=sqrt((muJ+muE)/(R1+R2)^3);

rI1=[-R1;0;0];
rI2=[R2;0;0];

statedot(1:3)=state(4:6); 
rI3=state(1:3);
r13=rI3-rI1;
r23=rI3-rI2;

statedot(4:6)=-muJ*r13/(norm(r13)^3)-muE*r23/(norm(r23)^3)-...
    2*cross([0;0;ws],state(4:6))-cross([0;0;ws],(cross([0;0;ws],rI3)));

end

function[statedot] = dVectJCI(t,state)
statedot=zeros(6,1); 
muJ=126.687e6; %km^3/s^2

statedot(1:3)=state(4:6); 
r = state(1:3);
statedot(4:6)=-muJ*r/(norm(r)^3);

end


function [r_out,v_out] = oeToSCI(a,e,i,Om,w,E,mu)
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

    