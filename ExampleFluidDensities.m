%% Compute Fluid Densities using Peng-Robinson EOS
% This script calls the computeFluidDensity function and plots the densities 
% of methane and CO2 as a function of pressure at a given temperature.
% Written by Filipe Borges in 2023

%% Preparing the Workspace
clear; 
close all; 
clc;

%% Set Graphic Properties
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultTextFontSize', 14); 
set(groot, 'DefaultTextFontName', 'Trebuchet MS');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultLegendFontSize', 14); 
set(groot, 'DefaultLegendFontName', 'Trebuchet MS');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultAxesFontSize', 12); 
set(groot, 'DefaultAxesFontName', 'Trebuchet MS');

%% Setup
% Define the range of pressures and the temperature for the analysis.
pressures = linspace(1e5, 3e7, 100); % Pressure range in Pascal
temperatureC = 35; % Temperature in Celsius

methaneDensities = zeros(size(pressures));
co2Densities = zeros(size(pressures));
methaneDensitiesIdealGas = zeros(size(pressures));

%% Compute Densities
% Call the computeFluidDensity function for each pressure value for each fluid.
for ii = 1:length(pressures)
    methaneDensities(ii) = computeFluidDensity(pressures(ii), temperatureC, "methane");
    co2Densities(ii) = computeFluidDensity(pressures(ii), temperatureC, "co2");
    methaneDensitiesIdealGas(ii) = computeMethaneDensityIdeal(pressures(ii), temperatureC);
end

%% Plotting
% Plot the densities as a function of pressure for both methane and CO2.

figure;
plot(pressures/1e5, methaneDensities, '-k', 'LineWidth', 2);
hold on;
plot(pressures/1e5, methaneDensitiesIdealGas, '-g', 'LineWidth', 2);
plot(pressures/1e5, co2Densities, '-r', 'LineWidth', 2);
hold off;
xlabel('Pressure (bar)');
ylabel('Density (kg/m$^3$)');
title(['Fluid Densities at ', num2str(temperatureC), ' $^{\circ}$C']);
legend('CH$_4$ - Peng-Robinson', 'CH$_4$ - Ideal Gas', 'CO$_2$ - Peng-Robinson', 'Location', 'East');
grid on;
