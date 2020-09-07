%% C4_parameter
%%Define the constants and parameters of the model

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fundamentals physics constants

R = 8.314;                          % Gas constant {J/K/mol]
RR = 1.987;                         % Gas constant [cal/K/mol]
Tk = 273.15;                        % Kelvin temperature constant [K]
c_e = 1.6021766208e-19;             % Elementary electric charge [C]
Kb = 1.38064852e-23;                % Boltzmann constant [J/K]
kb = 8.6173303e-5;                  % Boltzmann constant / c_e [eV/K]
Na = 6.022140857e23;                % Avogadro constant [/mol]
VV = 22.71;                         % Molar volume at STP [L/mol]


%% Zirconium material constants

MZr = 91.22;                        % Molar mass of zirconium [g/mol]
Mo = 16.0;                          % Molar mass of oxygen [g/mol]
MZi = MZr + 2*Mo;                   % Molar mass of zirconia [g/mol]
Mh = 1.0;                           % Molar mass of hydrogen [g/mol]
MNb = 93;                           % Molar mass of niobium [g/mol]
RhoZr = 6.509;                      % Density of zirconium [g/cm^3]
RhoZi = 5.68;                       % Density of zirconia [g/cm^3]
Vm = MZi /(RhoZi * Na);             % Average volume occupied by a molecule of zirconia [cm^3]
PBR = 1.55;                         % Pilling-Bedworth Ratio [] = molar volume of zirconia / molar volume of zirconium
T_transf = 950 + Tk;                % Temperature of phase transition alpha to beta [K]


%% Experimental parameters
%%Currently only used for hydrogen pickup

if strcmp(alloy, 'Zr05Nb')
    m_0 = 2387;                     % Initial mass of the sample [mg]
    Area = 16.5;                    % Surface area of the sample [cm^2]
    % Ch_nat = 14.37;                 % concentration of hydrogen in the natural material [wt ppm]
    Ch_nat = 14.12;                 % concentration of hydrogen in the natural material [wt ppm] (value given by the fit of the exp data)
elseif strcmp(alloy, 'Zr10Nb')
    m_0 = 1369;                     % Initial mass of the sample [mg]
    Area = 15.4;                    % Surface area of the sample [cm^2]
    % Ch_nat = 7.6;                   % concentration of hydrogen in the natural material [wt ppm]
    Ch_nat = 6.81;                  % concentration of hydrogen in the natural material [wt ppm] (value given by the fit of the exp data)
elseif strcmp(alloy, 'Zry4')        %% /!\ Not used because hydrogen pickup is usually neglected for high temperature corrosion of Zry-4
    m_0 = 2387;                     % Initial mass of the sample [mg]
    Area = 16.5;                    % Surface area of the sample [cm^2]
    % Ch_nat = 14.37;                 % concentration of hydrogen in the natural material [wt ppm]
    Ch_nat = 14.12;                 % concentration of hydrogen in the natural material [wt ppm] (value given by the fit of the exp data)
end

Nb_ox_state = 2;                    % Oxidation state of Nb to compensate the space charge []


%% Miscellaneous parameters

preoxide = 0e-4;                                    % Thickness of the pre-oxide layer [cm]

% Empirically determined oxide thickness at which the transition occured

if strcmp(mode, 'Monte Carlo') == 0     % If Monte Carlo transition study, the transition is defined in the Monte Carlo generation code
    if strcmp(alloy, 'Zr05Nb')
        transition = 1;             %% Current experimental data don't show any transition, so this value (1cm) is just too big to be reached
    elseif strcmp(alloy, 'Zr10Nb')
        transition = [3.1e-4 3.1e-4 3.1e-4 1];       % oxide thickness at the transition [cm] ; the "1" at the end is added artifically as a value too big to be reached
    else                            %% Notransition at high temperature
        transition = 1;
    end
end


%% LOCA temperature transient
%The exposure time is determined from the LOCA

% CINOG LOCA experiments
if strcmp(mode, 'LOCA')
    if strcmp(temperature, 'CINOG1')        % 20-700C @ 1C/s ; 700-1050C @ 0.1C/s
        T = [20:0.1:700, 700+0.01:0.01:1050]+Tk;
        exposure_time = '4180s';
    elseif strcmp(temperature, 'CINOG2')    % 20-1050C @ 1C/s
        T = (20:0.1:1050)+Tk;
        exposure_time = '1030s';
    elseif strcmp(temperature, 'CINOG3')    % 20-1200C @ 0.1C/s
        T = (20:0.01:1200)+Tk;
        exposure_time = '11800s';
    elseif strcmp(temperature, 'CINOG4')    % 20-1200C @ 1C/s
        T = (20:0.1:1200)+Tk;
        exposure_time = '1180s';
    elseif strcmp(temperature, 'CINOG5')    % 1050-1200C @ 0.1C/s
        T = (1050:0.01:1200)+Tk;
        exposure_time = '1500s';
    elseif strcmp(temperature, 'CINOG6')    % 1050-1200C @ 1C/s
        T = (1050:0.1:1200)+Tk;
        exposure_time = '150s';
    end
end

% ANL LOCA experiments
% The first 300 seconds is preoxidation at 360C, then
%The time spent at high temperature is used to determine the LOCA then the exposure time is redefined to include the entire time 

if strcmp(mode, 'LOCA') 
    if strcmp(temperature, '1200C') && strcmp(exposure_time, '300s')
%         T = [360*ones(1,3001) 1200*ones(1,3000)]+Tk;      % Use for simplified temperature transient
%         exposure_time = '600s';
        T = [360*ones(1,3000) 360:1.34:1200 1200*ones(1,12004)]+Tk;
        exposure_time = strcat(num2str((length(T)-1)/10), 's');
    elseif strcmp(temperature, '1200C') && strcmp(exposure_time, '600s')
%         T = [360*ones(1,3001) 1200*ones(1,6000)]+Tk;      % Use for simplified temperature transient
%         exposure_time = '900s';
        T = [360*ones(1,3000) 360:1.34:1200 1200*ones(1,12004)]+Tk;
        exposure_time = strcat(num2str((length(T)-1)/10), 's');
    elseif strcmp(temperature, '1200C') && strcmp(exposure_time, '1200s')
%         T = [360*ones(1,3001) 1200*ones(1,12000)]+Tk;     % Use for simplified temperature transient
%         exposure_time = '1500s';
        T = [360*ones(1,3000) 360:1.34:1200 1200*ones(1,12004)]+Tk;
        exposure_time = strcat(num2str((length(T)-1)/10), 's');
    elseif strcmp(temperature, '1100C') && strcmp(exposure_time, '1800s')
%         T = [360*ones(1,3001) 1100*ones(1,18000)]+Tk;     % Use for simplified temperature transient
%         exposure_time = '2100s';
        T = [360*ones(1,3000) 360:1.34:1200 1200*ones(1,12004)]+Tk;
        exposure_time = strcat(num2str((length(T)-1)/10), 's');
    end
end

% Leistikow LOCA experiments

if strcmp(mode, 'LOCA')
    if strcmp(temperature, 'Peak950C')
        T = [300:10:950, 950-2.5:-2.5:750]+Tk;
        exposure_time = '14.5s';
    elseif strcmp(temperature, 'Peak1000C')
        T = [300:10:1000, 1000-2.5:-2.5:750]+Tk;
        exposure_time = '17s';
        elseif strcmp(temperature, 'Peak1100C')
        T = [300:10:1100, 1100-2.5:-2.5:750]+Tk;
        exposure_time = '22s';
        elseif strcmp(temperature, 'Peak1200C')
        T = [300:10:1200, 1200-2.5:-2.5:750]+Tk;
        exposure_time = '27s';
        elseif strcmp(temperature, 'LeistikowLOCA')
        T = [300:10.833:950 950*ones(1,45) 950:-2.533:760 760*ones(1,60) 760:.517:1200 1200*ones(1,710) 1200:-8:800 800:-2.91:305]+Tk;
        exposure_time = strcat(num2str((length(T)-1)/10), 's');
    end
end


%% Time parameters

if exposure_time(end) == 'd'            % Time expressed in days
    t = str2double(exposure_time(1:end - 1)) * 86400;       % Duration of the simulation [s]
    N = 10 * str2double(exposure_time(1:end - 1));          % Number of time steps
elseif exposure_time(end) == 's'        % Time expressed in seconds
    t = str2double(exposure_time(1:end - 1));               % Duration of the simulation [s]
    N = 10 * str2double(exposure_time(1:end - 1));          % Number of time steps
end

dt = t / N;                         % time step duration [s]
time = linspace(0, t, N + 1);       % Vector containing the time values [s]

if exposure_time(end) == 'd'
    days = time / 86400;            % Vector containing the time values [days]
end


%% Isotherm temperature

if strcmp(mode, 'Isotherm') || strcmp(mode, 'Monte Carlo') || strcmp(mode, 'H optimization') || strcmp(mode, 'D HT optimization')
    if temperature(end) == 'C'      % Temperature expressed in degree Celsius
        T0 = str2double(temperature(1:end - 1)) + Tk;   % Isotherm temperature [K]
        T = T0 * ones(1, N + 1);        % Vector containing the temperature values at each time [K]
    elseif temperature(end) == 'K'  % Temperature expressed in Kelvin
        T0 = str2double(temperature(1:end - 1));        % Isotherm temperature [K]
        T = T0 * ones(1, N + 1);        % Vector containing the temperature values at each time [K]
    end
end


%% Linear temperature transient

if strcmp(mode, 'Linear')
    T1 = 700 + Tk;                  % Mid temperature [K]
    T2 = 1200 + Tk;                 % Final temperature [K]
    
    T = linspace(T1, T2, N + 1);    % Vector containing the temperature values at each time [K]
end

%% Space parameters
% The space mesh is divided into 2 different zones:
% - There is a fixed mesh far from the corroding surface
% - There is an adaptative mesh close to the corroding surface that moves with the interface
% The size of the adaptative mesh depends on the corrosion temperature
% A separate mesh is used in the oxide

e_Zr = 0.06;                                % Initial thickness of the cladding [cm]
if max(T) > T_transf
    e_fixed = 0.03;                         % thickness of the fixed mesh zone of the cladding [cm] (high temperature corrosion)
else
    e_fixed = 0.059;                        % thickness of the fixed mesh zone of the cladding [cm] (low temperature corrosion)
end
e_moving = e_Zr - e_fixed;                  % thickness of the adaptative mesh zone of the cladding [cm]
K_fixed = 100;                              % Number of spatial steps in the fixed mesh zone []
K_moving = 200;                             % Number of spatial steps in the adaptative mesh zone []
K_b = 100;                                  % Number of spatial steps in the beta phase []
K_a = K_moving - K_b;                       % Number of spatial steps in the alpha phase []
K = K_fixed + K_moving;                     % Total number of spatial steps []
K_ox = 100;                                 % Number of spatial steps in the oxide []


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Constants related to C4 model

nu = 1e13;          % Frequency of jump migration [Hz]
a = 5e-8;           % length of the jump, lattice parameter [cm]

Zv = 2;             % Number of charge of vacancies []
Ze = -1;            % Number of charges of electrons []
Zh = 1;             % Number of charges of hydrogen []

%% Migration energies
%Migration energies depend on the alloy, and on the model used (see comments in file C4_launch)
%Emv is the migration energy of vacancies [eV]
%Eme is the migration energy of electrons [eV]
%Emh is the migration energy of hydrogen [eV]

%% Positive space charge: Emv > Eme

if strcmp(alloy, 'Zr05Nb')
    if strcmp(model, 'C4')
        Emv = 1.53 * ones(1, N + 1);
        Eme = 1.35 * ones(1, N + 1);
    elseif strcmp(model, 'C4-O')
        Emv = 1.53 * ones(1, N + 1);
        Eme = 1.35 * ones(1, N + 1);
    elseif strcmp(model, 'C4-H')
        Emv = 1.53 * ones(1, N + 1);
        Eme = 1.32 * ones(1, N + 1);
        Emh = 1.5 * ones(1, N+ 1);                          % Constant hydrogen migration energy
%         load('input_file/Emh_Zr05Nb_C4-H_360C_240d');       % Optimized hydrogen migration energy to fit H concentration at operating temperature
    elseif strcmp(model, 'C4-OH')
        Emv = 1.53 * ones(1, N + 1);
        Eme = 1.35 * ones(1, N + 1);
        Emh = 1.7 * ones(1, N + 1);                         % Constant hydrogen migration energy
%         load('input_file/Emh_Zr05Nb_C4-OH_360C_240d');      % Optimized hydrogen migration energy to fit H concentration at operating temperature (TO BE UPDATED)
    end
elseif strcmp(alloy, 'Zr10Nb')
    if strcmp(model, 'C4')
        Emv = 1.43 * ones(1, N + 1);
        Eme = 1.35 * ones(1, N + 1);
    elseif strcmp(model, 'C4-O')
        Emv = 1.43 * ones(1, N + 1);
        Eme = 1.35 * ones(1, N + 1);
    elseif strcmp(model, 'C4-H')
        Emv = 1.43 * ones(1, N + 1);
        Eme = 1.35 * ones(1, N + 1);
        Emh = 1.65 * ones(1, N + 1);                        % Constant hydrogen migration energy
%         load('input_file/Emh_Zr10Nb_C4-H_360C_240d');       % Optimized hydrogen migration energy to fit H concentration at operating temperature
    elseif strcmp(model, 'C4-OH')
        Emv = 1.47 * ones(1, N + 1);
        Eme = 1.3 * ones(1, N + 1);
        Emh = 1.65 * ones(1, N + 1);                        % Constant hydrogen migration energy
%         load('input_file/Emh_Zr10Nb_C4-OH_360C_240d');      % Optimized hydrogen migration energy to fit H concentration at operating temperature (TO BE UPDATED)
    end
elseif max(T) >= T_transf
    if strcmp(model, 'C4') || strcmp(model, 'C4-H')
        error('Oxygen dissolution is required for high-temperature models');
    elseif strcmp(model, 'C4-O')
        Emv = 1.35 * ones(1, N + 1);
        Eme = 1.30 * ones(1, N + 1);
    elseif strcmp(model, 'C4-OH')       % NOT UP TO DATE
        Emv = 1.57 * ones(1, N + 1);
        Eme = 1.0 * ones(1, N + 1);
        Emh = 1.6 * ones(1, N + 1);
    end
elseif strcmp(alloy, 'Zry4') && strcmp(temperature, '950C')       % NOT UP TO DATE
    if strcmp(model, 'C4') || strcmp(model, 'C4-H')
        error('Oxygen dissolution is required for high-temperature models');
    elseif strcmp(model, 'C4-O')
        Emv = 1.47 * ones(1, N + 1);
        Eme = 1.3 * ones(1, N + 1);
    elseif strcmp(model, 'C4-OH')
        Emv = 1.57 * ones(1, N + 1);
        Eme = 1.0 * ones(1, N + 1);
        Emh = 1.6 * ones(1, N + 1);
    end
end

if strcmp(model, 'C4') || strcmp(model, 'C4-O')
    Emh = 0 * ones(1, N + 1);
end


%% Concentration at each interface
%Only concentrations independant of the temperature are defined here
%Temperature-dependent concentration are defined in the C4_oxidation code

% Oxide-water interface concentrations

Cv_ox_w = 1.0e17;               % Concentration of vacancies in the oxide at the oxide/water interface [at/cm^3]
Ce_ox_w = 2.0e17;               % Concentration of electrons in the oxide at the oxide/water interface [at/cm^3]
Ch_ox_w = 3e21 * ones(1, N+1);  % Concentration of hydrogen in the oxide at the oxide/water interface [at/cm^3]
Ch_ox_w_ppm = Ch_ox_w * Mh / (RhoZi * Na) * 1e6;        % Concentration of hydrogen in the oxide at the oxide/water interface [wt ppm]
% Ch_ox_w = 1.8e21;               % Concentration of hydrogen in the oxide at the oxide/water interface [at/cm^3] (C4+ model)

% Oxide-metal interface concentrations

xo_stoichio = 2/3;
Co_stoichio = xo_stoichio / (1 - xo_stoichio) * RhoZi * Na / MZi;               % Concentration of oxygen in the oxide at the oxide/alpha-phase interface [at/cm^3]

Ch_ox_a = Ch_nat * 1e-6 * RhoZi * Na / Mh;               % Concentration of hydrogen in the oxide at the oxide/alpha-phase interface [at/cm^3]     %% Use RhoZi instead of RhoZr because in the oxide

% Metal-oxide interface concentrations

if strcmp(model, 'C4') || strcmp(model, 'C4-H')
    Co_a_ox = 0;                                            % Concentration of oxygen in the alpha-phase at the oxide/alpha-phase interface [at/cm^3]
    xo_a_ox = 1 / (1 + RhoZr * Na / (MZr .* Co_a_ox));      % Fraction of oxygen in the alpha-phase at the oxide/alpha-phase interface []
    Co_nat = 0;                                             % Concentration of oxygen in the natural material [at/cm^3]
    xo_nat = 0;                                             % Fraction of oxygen in the natural material []
elseif strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')
    xo_nat = 0.0074;                                            % Oxygen atomic fraction in the natural material []
%     xo_nat = 0.013;                                             % Oxygen atomic fraction in the natural material (from EPMA) []
    Co_nat = xo_nat / (1 - xo_nat) * RhoZr * Na / MZr;          % Concentration of oxygen in the natural material [at/cm^3] (4.9e20)
    xo_duct_b = 0.022;                                          % Fraction of oxygen at which the beta phase becomes ductile [] (From Brachet)
    Co_duct_b = xo_duct_b / (1 - xo_duct_b) * RhoZr * Na / MZr; % Concentration of oxygen at which the beta phase becomes ductile [at/cm^3]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save all these variables in the input file
save(strcat('input_file/',file_name), '-regexp', '^(?!(progress)$).');