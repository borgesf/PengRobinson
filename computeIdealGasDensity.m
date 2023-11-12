function rho = computeIdealGasDensity(P, T_Celsius, gasType)
    %% Computes gas density using the Ideal Gas Law
    % Based on the Ideal Gas equation: PV = nRT
    % Written by Filipe Borges in 2023
    % INPUT:
    %   - P: Pressure in Pascal
    %   - T_Celsius: Temperature in Celsius
    %   - gasType: Either 'ch4' (methane) or 'co2'
    % OUTPUT:
    %   - rho: Gas density in kg/m^3

    arguments
        P (1,1) double {mustBePositive}
        T_Celsius (1,1) double
        gasType (1,1) string {mustBeMember(gasType, ["ch4", "co2"])}
    end

    R = 8.314; % Universal gas constant in J/(mol*K)
    T = T_Celsius + 273.15; % Convert temperature to Kelvin

    Mw = getGasMolarWeight(gasType); % Molar weight in kg/mol

    %% Calculate density using rho = P*Mw / (R*T)
    rho = (P * Mw) / (R * T); % density in kg/m^3
end

function Mw = getGasMolarWeight(gasType)
    %% Compute the molar weight of the specified gas
    % Written by Filipe Borges in 2023
    % INPUT:
    %   - gasType: Either 'CH4' or 'CO2'
    % OUTPUT:
    %   - Mw: Molar weight of the gas in kg/mol

    switch gasType
        case "ch4"
            Mw = 0.01604; % kg/mol for methane
        case "co2"
            Mw = 0.04401; % kg/mol for CO2
        otherwise
            error("Invalid gas type. Choose either 'ch4' or 'co2'.");
    end
end
