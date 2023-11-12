function rho = computeDensityPengRobinson(P, T_Celsius, fluid)
    %% Computes the fluid density using Peng-Robinson EOS
    % Based on the Peng-Robinson equation of state, as described in: 
    % "Peng, D. Y., & Robinson, D. B. (1976). A New Two-Constant Equation of State. 
    % Industrial & Engineering Chemistry Fundamentals, 15(1), 59–64."
    % Written by Filipe Borges in 2023
    % INPUT:
    %   - P: Pressure in Pascal
    %   - T_Celsius: Temperature in Celsius
    %   - fluid: Either 'ch4' or 'co2'
    % OUTPUT:
    %   - rho: Fluid density in kg/m^3

    arguments
        P (1,1) double {mustBePositive}
        T_Celsius (1,1) double
        fluid (1,1) string {mustBeMember(fluid, ["ch4", "co2"])}
    end

    R = 8.314; % Universal gas constant in J/(mol*K)
    T = T_Celsius + 273.15; % Convert temperature to Kelvin

    switch fluid
        case "ch4"
            Tc = 190.56; % Critical temperature (K)
            Pc = 4.59e6; % Critical pressure (Pa)
            omega = 0.011;
        case "co2"
            Tc = 304.13; % Critical temperature (K)
            Pc = 7.38e6; % Critical pressure (Pa)
            omega = 0.225;
        otherwise
            error("Invalid fluid type. Choose either 'ch4' or 'co2'.");
    end

    a = 0.45724 * (R^2 * Tc^2) / Pc;
    b = 0.07780 * R * Tc / Pc;

    alpha = (1 + (0.37464 + 1.54226 * omega - 0.26992 * omega^2) * (1 - sqrt(T / Tc)))^2;
    A = alpha * a * P / (R^2 * T^2);
    B = b * P / (R * T);

    %% Cubic equation: Z^3 + c2*Z^2 + c1*Z + c0 = 0
    c2 = -(1 - B);
    c1 = A - 3 * B^2 - 2 * B;
    c0 = -(A * B - B^2 - B^3);

    % Solve for Z
    Z_roots = roots([1, c2, c1, c0]);
    Z_roots_real = Z_roots(imag(Z_roots) == 0);

    % Choose the largest Z for gas phase and smallest for liquid
    Z = max(Z_roots_real);

    %% Calculate molar volume using Z = PV/RT
    Vm = Z * R * T / P; % molar volume in m^3/mol
    Mw = getMolarWeight(fluid); % Molar weight in kg/mol
    rho = Mw / Vm; % density in kg/m^3
end

function Mw = getMolarWeight(fluid)
    %% Compute the molar weight based on the fluid type
    % Written by FIBO in 2023
    % INPUT:
    %   - fluid: Either 'methane' or 'co2'
    % OUTPUT:
    %   - Mw: Molar weight in kg/mol

    switch fluid
        case "ch4"
            Mw = 0.01604; % kg/mol
        case "co2"
            Mw = 0.04401; % kg/mol
    end
end
