%function drag = bEuropaAtmosphericModel(height,velocity)
height = 1000; 
velocity = 1.312e3;
g = 1.315; %m/s^2
gamma = 1.4; %O2 atmosphere, diatomic value
Ru = 8314.5; 
MW = 32; %g/cm^3

numDens = 1000e8; %number density of oxygen on surface, particles/cm^3
mols = 6.02e23; %avogadro's number

surfaceDensity = numDens*MW/mols*1000; %surface density in kg/m^3
R = Ru/MW; 
T0 = 100; %K

Cd = 1.5; 
radius = 2.65/2; %m
A = pi*radius^2; 

density = surfaceDensity*(1-(gamma-1)*g*height/(gamma*R*T0))^(1/(gamma-1)); 

drag = Cd*density*velocity^2*A/2; 
%fprintf('Drag = %f N\n', drag)
%end
