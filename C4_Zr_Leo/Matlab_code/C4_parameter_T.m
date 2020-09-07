%% C4_parameter_T
%%Function computing the temperature-dependent parameters

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated: 09/08/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute the temperature-dependent parameters

if strcmp(temperature, '360C')
    substoichio = 0.0076;                               % From discussion with Corvolan
    xo_ox_a = 2/3 * (1 - substoichio);
else
    xo_ox_a = -1.10606E-05 * T(n) + 0.667118123;        % From phase diagram (Ma 2008)
end
Co_ox_a = xo_ox_a * 3 * RhoZi * Na / MZi;               % Concentration of oxygen in the oxide at the oxide/alpha-phase interface [at/cm^3]
substoichio = 1 - xo_ox_a * 3 / 2;
Cv_ox_a = Co_stoichio - Co_ox_a;                        % Concentration of vacancies in the oxide at the oxide/alpha-phase interface [at/cm^3]
Ce_ox_a = 2 * Cv_ox_a;                                  % Concentration of electrons in the oxide at the oxide/alpha-phase interface [at/cm^3]

% Set up the oxygen concentration in the metal at the oxide/alpha interface and at the alpha/beta interface
%Equations come from Abriata, 1986
if strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')
    if T(n) >= 273.15 && T(n) <= 473.15           % Artificially added (not present in the Zr-O phase diagram)
        xo_a_ox = 28.6 * 1e-2;
    elseif T(n) >= 473.15 && T(n) <= 1478.15
        xo_a_ox = (28.6 + exp(-6748/T(n) + 4.748)) * 1e-2;
    elseif T(n) > 1478.15 && T(n) <= 1798.15
        xo_a_ox = (28.6 + exp(-6301/T(n) + 4.460)) * 1e-2;
    elseif T(n) > 1798.15 && T(n) <= 2338.15
        xo_a_ox = (28.6 + exp(-7012/T(n) + 8.434 - 3.521e-3*T(n) + 0.851e-6*T(n)^2)) * 1e-2;
    end
    Co_a_ox = xo_a_ox / (1 - xo_a_ox) * RhoZr * Na / MZr;               % Concentration of oxygen in the alpha-phase at the oxide/alpha-phase interface [at/cm^3]
    xo_a_b = (45.86e-3*(T(n)-1136.15) - 44.77e-6*(T(n)-1136.15)^2 + 17.40e-9*(T(n)-1136.15)^3) * 1e-2;
    Co_a_b = xo_a_b * RhoZr * Na / (1 - xo_a_b) / MZr;                  % Concentration of oxygen in the alpha-phase at the alpha-/beta-phase interface [at/cm^3]
    xo_b_a = (9.59e-3*(T(n)-1136.15) + 4.72e-6*(T(n)-1136.15)^2 - 4.35e-9*(T(n)-1136.15)^3) * 1e-2;
    Co_b_a = xo_b_a * RhoZr * Na / (1 - xo_b_a) / MZr;                  % Concentration of oxygen in the beta-phase at the alpha-/beta-phase interface [at/cm^3]
end