% AA279B - Class Project
% From Earth to the Water on Europa
% Annex 1: Low Orbit Definition 
% Code by:
% Christine Hamilton and
% Juan Martinez Castellanos
clc; clear all; close all;

G = 6.67430*10^-11; % [m3/kg1s2]
% Earth
Me = 5.9724*10^24; % [kg]
Re = 6371000; % [m]
LEO = 2000000; % [m]
% Jupiter
Mj = 1898.19*10^24; % [kg]
Rj = 69911000; % [m]
% Europa
Meu = 48*10^21; % [kg]
Reu = 3122000/2; % [m]

acc = G*Me/(LEO+Re)^2;
LJO = sqrt(G*Mj/acc)-Rj;
LEuO = sqrt(G*Meu/acc)-Reu;
fprintf("Low Orbit Acceleration = %f km/s2\n",acc/1000);
fprintf("Low Earth Orbit = %f km\n",LEO/1000);
fprintf("Low Jupiter Orbit = %f km\n",LJO/1000);
fprintf("Low Europa Orbit = %f km\n",LEuO/1000);
