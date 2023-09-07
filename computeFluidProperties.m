function [rho, c] = computeFluidProperties(P, temperatureCelsius, fluidType)
    %% Computes the fluid density and speed of sound using Peng-Robinson EOS
    % Based on the Peng-Robinson equation of state, as described in: 
    % "Peng, D. Y., & Robinson, D. B. (1976). A New Two-Constant Equation of State. 
    % Industrial & Engineering Chemistry Fundamentals, 15(1), 59–64."
    % Written by Filipe Borges in 2023
    % INPUT:
    %   - P: Pressure in Pascal
    %   - temperatureCelsius: Temperature in Celsius
    %   - fluidType: Either 'CH4' or 'CO2'
    % OUTPUT:
    %   - rho: Fluid density in kg/m^3
    %   - c: Speed of sound in the fluid (m/s)

    arguments
        P (1,1) double {mustBePositive}
        temperatureCelsius (1,1) double
        fluidType (1,1) string {mustBeMember(fluidType, ["CH4", "CO2"])}
    end

    %% Constants
    universalGasConstant = 8.314; % Universal gas constant in J/(mol*K)
    temperatureKelvin = temperatureCelsius + 273.15; % Convert temperature to Kelvin

    %% Determine fluid constants based on type
    switch fluidType
        case "CH4"
            criticalTemperature = 190.56; % Critical temperature (K)
            criticalPressure = 4.59e6; % Critical pressure (Pa)
            acentricFactor = 0.011;
        case "CO2"
            criticalTemperature = 304.13; % Critical temperature (K)
            criticalPressure = 7.38e6; % Critical pressure (Pa)
            acentricFactor = 0.225;
        otherwise
            error("Invalid fluid type. Choose either 'CH4' or 'CO2'.");
    end

    attractionParam = 0.45724 * (universalGasConstant^2 * criticalTemperature^2) / criticalPressure;
    excludedVolume = 0.07780 * universalGasConstant * criticalTemperature / criticalPressure;

    alpha = (1 + (0.37464 + 1.54226 * acentricFactor - 0.26992 * acentricFactor^2) * (1 - sqrt(temperatureKelvin / criticalTemperature)))^2;
    coefA = alpha * attractionParam * P / (universalGasConstant^2 * temperatureKelvin^2);
    coefB = excludedVolume * P / (universalGasConstant * temperatureKelvin);

    %% Cubic equation coefficients
    coefZ2 = -(1 - coefB);
    coefZ1 = coefA - 3 * coefB^2 - 2 * coefB;
    coefZ0 = -(coefA * coefB - coefB^2 - coefB^3);

    % Solve for compressibility factor Z
    Z_roots = roots([1, coefZ2, coefZ1, coefZ0]);
    Z_roots_real = Z_roots(imag(Z_roots) == 0);

    % Choose the largest Z for gas phase
    compressibilityFactor = max(Z_roots_real);

    %% Calculate molar volume and density
    molarVolume = compressibilityFactor * universalGasConstant * temperatureKelvin / P;
    molarWeight = getMolarWeight(fluidType); 
    rho = molarWeight / molarVolume; 

    %% Speed of sound calculations (using derivatives of Peng-Robinson EOS)
    dPdV = -P/molarVolume^2 + (2*attractionParam/molarVolume^3)*(1-coefB/molarVolume) ...
           - attractionParam*coefB/molarVolume^4 ...
           - (P*coefB+attractionParam*coefB^2)/molarVolume^2;
    c = sqrt(-dPdV * molarWeight / rho^2);
end

function molarWeight = getMolarWeight(fluidType)
    %% Compute the molar weight based on the fluid type
    % Written by Filipe Borges in 2023
    % INPUT:
    %   - fluidType: Either 'CH4' or 'CO2'
    % OUTPUT:
    %   - molarWeight: Molar weight in kg/mol

    switch fluidType
        case "CH4"
            molarWeight = 0.01604; % kg/mol
        case "CO2"
            molarWeight = 0.04401; % kg/mol
    end
end
