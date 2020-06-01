function [E]=meanToEcc(M,e,tol)
%M and tol in radians
if e>0.95
    E=M;
else
    E=pi;
end
d=-(E-e*sin(E)-M)/(1-e*cos(E));
while abs(d)>tol
    d=-(E-e*sin(E)-M)/(1-e*cos(E));
    E=E+d;
end
