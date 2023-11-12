%% Compute Fluid Densities using Peng-Robinson EOS
% This script calls the computeFluidDensity function and plots the densities 
% of methane and CO2 as a function of pressure at a given temperature.
% Written by Filipe Borges in 2023

%% Setup
% Define the range of pressures and the temperature for the analysis.
pressures = linspace(1e5, 3e7, 500); % Pressure range in Pascal
temperatureC = 35; % Temperature in Celsius

methaneDensitiesPR = zeros(size(pressures));
co2DensitiesPR = zeros(size(pressures));
methaneDensitiesIdealGas = zeros(size(pressures));
co2DensitiesIdealGas = zeros(size(pressures));

%% Compute Densities
% Call the computeFluidDensity function for each pressure value for each fluid.
for ii = 1:length(pressures)
    methaneDensitiesPR(ii) = computeDensityPengRobinson(pressures(ii), temperatureC, "ch4");
    co2DensitiesPR(ii) = computeDensityPengRobinson(pressures(ii), temperatureC, "co2");
    methaneDensitiesIdealGas(ii) = computeDensityIdealGas(pressures(ii), temperatureC,"ch4");
    co2DensitiesIdealGas(ii) = computeDensityIdealGas(pressures(ii), temperatureC,"co2");
end

%% Plotting
% Plot the densities as a function of pressure for both methane and CO2.

figure;
set(gcf, 'Position', [100, 100, 1400, 600]);
plot(pressures/1e5, methaneDensitiesPR, '-k', 'LineWidth', 2,'DisplayName','CH$_4$ - Peng-Robinson');
hold on;
box on;
plot(pressures/1e5, methaneDensitiesIdealGas, '-g', 'LineWidth', 2,'DisplayName','CH$_4$ - Ideal gas');
plot(pressures/1e5, co2DensitiesPR, '-r', 'LineWidth', 2,'DisplayName','CO$_2$ - Peng-Robinson');
plot(pressures/1e5, co2DensitiesIdealGas, '-b', 'LineWidth', 2,'DisplayName','CO$_2$ - Ideal gas');
hold off;
xlabel('Pressure (bar)');
ylabel('Density (kg/m$^3$)');
title(['Fluid Densities at ', num2str(temperatureC), ' $^{\circ}$C']);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
legend('Location', 'northwest');
grid on;
