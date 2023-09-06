function rho = computeMethaneDensityIdeal(P, T_Celsius)
    %% Computes methane density using the Ideal Gas Law
    % Based on the Ideal Gas equation: PV = nRT
    % Written by Filipe Borges in 2023
    % INPUT:
    %   - P: Pressure in Pascal
    %   - T_Celsius: Temperature in Celsius
    % OUTPUT:
    %   - rho: Methane density in kg/m^3

    arguments
        P (1,1) double {mustBePositive}
        T_Celsius (1,1) double
    end

    R = 8.314; % Universal gas constant in J/(mol*K)
    T = T_Celsius + 273.15; % Convert temperature to Kelvin

    Mw = getMethaneMolarWeight(); % Molar weight in kg/mol

    %% Calculate density using rho = P*Mw / (R*T)
    rho = (P * Mw) / (R * T); % density in kg/m^3
end

function Mw = getMethaneMolarWeight()
    %% Compute the molar weight of methane
    % Written by FIBO in 2023
    % OUTPUT:
    %   - Mw: Molar weight of methane in kg/mol

    Mw = 0.01604; % kg/mol for methane
end
