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

ch4DensityPengRob = zeros(size(pressures));
ch4DensityIdealGas = zeros(size(pressures));
ch4DensityvanDerWaals = zeros(size(pressures));
ch4SpeedPengRob = zeros(size(pressures));
ch4SpeedvanDerWaals = zeros(size(pressures));

co2DensityPengRob = zeros(size(pressures));
co2DensityvanDerWaals = zeros(size(pressures));
co2SpeedPengRob = zeros(size(pressures));
co2SpeedvanDerWaals = zeros(size(pressures));


%% Compute Densities
% Call the computeFluidDensity function for each pressure value for each fluid.
for ii = 1:length(pressures)
    [ch4DensityPengRob(ii), ch4SpeedPengRob(ii)] = computeFluidProperties(pressures(ii), temperatureC, "CH4");
    ch4DensityIdealGas(ii) = computeMethaneDensityIdeal(pressures(ii), temperatureC);
    [ch4DensityvanDerWaals(ii), ch4SpeedvanDerWaals(ii)] = computeGasVanDerWaals(pressures(ii), temperatureC,'CH4');

    [co2DensityPengRob(ii), co2SpeedPengRob(ii)] = computeFluidProperties(pressures(ii), temperatureC, "CO2");
    [co2DensityvanDerWaals(ii), co2SpeedvanDerWaals(ii)] = computeGasVanDerWaals(pressures(ii), temperatureC,'CO2');
end

%% Plotting
% Plot the densities as a function of pressure for both methane and CO2.

figure;
set(gcf, 'Position', [0, 100, 1200, 500]);
subplot(1, 2, 1);
plot(pressures/1e5, ch4DensityPengRob, '-k', 'LineWidth', 2);
hold on;
box on;
plot(pressures/1e5, ch4DensityvanDerWaals, '-c', 'LineWidth', 2);
plot(pressures/1e5, ch4DensityIdealGas, '-g', 'LineWidth', 2);
plot(pressures/1e5, co2DensityPengRob, '-r', 'LineWidth', 2);
plot(pressures/1e5, co2DensityvanDerWaals, '-b', 'LineWidth', 2);
hold off;
xlabel('Pressure (bar)');
ylabel('Density (kg/m$^3$)');
title(['Fluid Density at ', num2str(temperatureC), ' $^{\circ}$C']);
legend('CH$_4$ - Peng-Robinson', 'CH$_4$ - van Der Waals', 'CH$_4$ - Ideal Gas',...
       'CO$_2$ - Peng-Robinson','CO$_2$ - van Der Waals', 'Location', 'East');
grid on;

subplot(1, 2, 2);
plot(pressures/1e5, ch4SpeedPengRob, '-k', 'LineWidth', 2);
hold on;
plot(pressures/1e5, ch4SpeedvanDerWaals, '-c', 'LineWidth', 2);
box on;
%plot(pressures/1e5, ch4DensityIdealGas, '-g', 'LineWidth', 2);
plot(pressures/1e5, co2SpeedPengRob, '-r', 'LineWidth', 2);
plot(pressures/1e5, co2SpeedvanDerWaals, '-b', 'LineWidth', 2);
hold off;
xlabel('Pressure (bar)');
ylabel('Speed of Sound (m/s)');
title(['Speed of Sound at ', num2str(temperatureC), ' $^{\circ}$C']);
%legend('CH$_4$ - Peng-Robinson', 'CH$_4$ - van Der Waals', 'CH$_4$ - Ideal Gas',...
%       'CO$_2$ - Peng-Robinson','CO$_2$ - van Der Waals', 'Location', 'East');
grid on;
