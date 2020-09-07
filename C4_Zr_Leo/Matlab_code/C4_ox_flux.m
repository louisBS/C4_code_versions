%% C4_ox_flux
%%Function computing the constant used to calculate the different fluxes (vacancy, electron and hydrogen) in the oxide

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated: 09/08/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute the mobility of the different species

mu_v = 4 * a^2 * nu * Zv / (kb * T(n)) * exp(-Emv(n) / (kb * T(n)));      % Mobility of vacancies [cm^2/V/s]
mu_e = 4 * a^2 * nu * Ze / (kb * T(n)) * exp(-Eme(n) / (kb * T(n)));      % Mobility of electrons [cm^2/V/s]

D_v = 4 * a^2 * nu * exp(-Emv(n) / (kb * T(n)));      % Vacancy diffusion coefficient [cm^2/s]
D_e = 4 * a^2 * nu * exp(-Eme(n) / (kb * T(n)));      % Electron diffusion coefficient [cm^2/s]

%% Define the terms of the second order equation that is solved.
%%The form of the equation is A*x^2 + B*x + C = 0.
%Different cases, depending on the model choice (hydrogen transport or not)

if strcmp(model, 'C4') || strcmp(model, 'C4-O')         % Models without hydrogen transport
    A = mu_e * Ce_ox_w - 2 * mu_v * Cv_ox_a;
    B = mu_e * Ce_ox_w - mu_e * Ce_ox_a;
    C = 2 * mu_v * Cv_ox_w - mu_e * Ce_ox_a;
elseif strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')  % Models with hydrogen transport
    mu_h = 4 * a^2 * nu * Zh / (kb * T(n)) * exp(-Emh(n) / (kb * T(n)));  % Mobility of hydrogen [cm^2/V/s]
    D_h(n) = 4 * a^2 * nu * exp(-Emh(n) / (kb * T(n)));      % Hydrogen diffusion coefficient [cm^2/s]
    
    A = mu_e * Ce_ox_w - 2 * mu_v * Cv_ox_a - mu_h * Ch_ox_a;
    B = mu_e * Ce_ox_w - mu_e * Ce_ox_a + mu_h * Ch_ox_w(n) - mu_h * Ch_ox_a;
    C = 2 * mu_v * Cv_ox_w - mu_e * Ce_ox_a + mu_h * Ch_ox_w(n);
end

delta = B^2 - 4 * A * C;                % Determinant of the equation
eta = (-B - sqrt(delta)) / (2 * A);     % Solution of the equation

potential = kb * T(n) * log(eta);      % Potential [V]

gammav = mu_v * potential * (Cv_ox_w - Cv_ox_a * eta^2)/(1 - eta^2);    % Constant used for the calculation of the vacancy flux
gammae = mu_e * potential * (Ce_ox_w * eta - Ce_ox_a)/(eta - 1);        % Constant used for the calculation of the electron flux

if strcmp(model, 'C4') || strcmp(model, 'C4-O')         % Models without hydrogen transport
    Ch_ox_a = 0;        % Set the concentration of hydrogen at the oxide-metal interface to 0
    gammah = 0;         % Set the hydrogen flux to 0
elseif strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')  % Model with hydrogen transport
    gammah = mu_h * potential * (Ch_ox_a * eta - Ch_ox_w(n))/(1 - eta);    % Constant used for the calculation of the hydrogen flux
end
