%% C4_isotherm_spinel
% Calculates oxide thickness for isothermal corrosion of Cr and Mn rich
% steel.
% One layer of oxide is modelled : MnCr2O4 spinel. It is modelled as a
% p-type oxide with diffusion of cations (vacancies) and outward growth.

%%Author : Louis Bailly-Salins
%%Email : baillysalins@wisc.edu

%%Last updated : 30/07/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear all;
close all;
%% Fundamentals physics constants

R = 8.314;                          % Gas constant {J/K/mol]
RR = 1.987;                         % Gas constant [cal/K/mol]
Tk = 273.15;                        % Kelvin temperature constant [K]
c_e = 1.6021766208e-19;             % Elementary electric charge [C]
Kb = 1.38064852e-23;                % Boltzmann constant [J/K]
kb = 8.6173303e-5;                  % Boltzmann constant / c_e [eV/K]
Na = 6.022140857e23;                % Avogadro constant [/mol]
VV = 22.71;                         % Molar volume at STP [L/mol]

%% Material parameters

%Composition of the steel (wt% / at%) :
% Fe : 67.82 / 65.48
% Cr : 20.25 / 21.00
% Mn : 8.5 / 8.34
% Ni : 2.13 / 1.96
% C  : 0.55 / 2.47
% Mo : 0.5 / 0.28
% Si : 0.25 / 0.48

x_steel_Fe = 0.6548;                % Native atomic fraction of Fe in 21-2N steel [at %]
x_steel_Cr = 0.2100;                % Native atomic fraction of Cr in 21-2N steel [at %]
x_steel_Mn = 0.0834;                % Native atomic fraction of Mn in 21-2N steel [at %]
x_steel_Ni = 0.0196;                % Native atomic fraction of Ni in 21-2N steel [at %]
x_steel_C = 0.0247;                 % Native atomic fraction of C in 21-2N steel [at %]
x_steel_Mo = 0.0028;                % Native atomic fraction of Mo in 21-2N steel [at %]
x_steel_Si = 0.0048;                % Native atomic fraction of Si in 21-2N steel [at %]

MFe = 55.845;                       % Molar mass of Fe [g/mol]
MMn = 54.938044;                    % Molar mass of Mn [g/mol]
MCr = 51.9961;                      % Molar mass of Cr [g/mol]
MNi = 58.6934;                      % Molar mass of Ni [g/mol]
MC = 12.0107;                       % Molar mass of C [g/mol]
MMo = 95.95;                        % Molar mass of Mo [g/mol]
MSi = 28.0855;                      % Molar mass of Si [g/mol]
Msteel = x_steel_Fe*MFe + x_steel_Cr*MCr + x_steel_Mn*MMn + x_steel_Ni*MNi + x_steel_C*MC + x_steel_Mo*MMo + x_steel_Si*MSi;

MO = 15.999;                               % Molar mass of oxygen [g/mol]
Mspinel = MMn + 2*MCr + 4*MO;              % Molar mass of zirconia [g/mol]

Rhosteel = 7.67;                           % Density of 21-2N steel [g/cm^3]
Rhospinel = 4.93;                          % Density of zirconia [g/cm^3]

Vspinel = 1e21 * Mspinel/(Rhospinel * Na); % Average volume occupied by a molecule of MnCr2O4 [nm^3]
Cspinel = 1/Vspinel;                       % Concentration of MnCr2O4 molecules [molecules/nm^3]   

Csteel = 1e-21 * Rhosteel*Na/Msteel ;      % Concentration of metal atoms in steel 21-2N [at/nm^3]

%% C4 model parameters

%nu = 1e13;               % Frequency of jump migration [Hz]
nu = 3600 * 1e13;         % Frequency of jump migration [/hr]
a = 0.5;                  % length of the jump, lattice parameter [nm]

ZVMn = -2;                % Number of charge of Mn vacancies []
Zh = 1;                   % Number of charges of holes []

%% Migration energies 

EmVMn = 1.82;        % Migration energy of Mn vacancies [eV]
Emh = 1.7;          % Migration energy of electron holes [eV]

%% Simulation parameters : Temperature and time 

T = 700 + 273.15 ;           % Temperature [K]
%t = 3600 * 150              % Duration of experiment [s]
t = 500;                      % Duration of experiment [hr]

%% Interfaces concentrations

% Oxide-gas interface concentrations
C_VMn_ox_g = 1e-4 * Cspinel ;          % Concentration of Mn vacancies in the oxide at the oxide/gas interface [at/cm^3]   
C_VCr_ox_g = 1e-4 * 2 * Cspinel;
C_h_ox_g = 2*C_VMn_ox_g + 3*C_VCr_ox_g ; % Concentration of holes in the oxide at the oxide/gas interface [at/cm^3] 

% Oxide-metal interface concentrations (negligible)
C_VMn_ox_m = 1e-8;                       % Concentration of Mn vacancies in the oxide at the oxide/metal interface [at/cm^3]
C_h_ox_m = 2e-8;                         % Concentration of holes in the oxide at the oxide/metal interface [at/cm^3]

%% Mobilities

mu_VMn = 4*a^2*nu*ZVMn/(kb*T)*exp(-EmVMn/(kb*T));  % Mobilty of Mn vacancies [nm^2/V/hr]
mu_h = 4*a^2*nu*Zh/(kb*T)*exp(-Emh/(kb*T));        % Mobilty of holes [nm^2/V/hr]

%% Solving for eta

A = 8*mu_VMn*C_VMn_ox_g - mu_h*C_h_ox_m;
B = mu_h*(C_h_ox_g-C_h_ox_m);
C = mu_h*C_h_ox_g - 8*mu_VMn*C_VMn_ox_m;

eta = (-B-sqrt(B^2-4*A*C))/(2*A);

%% Computing gamma and parabolic growth rate constant kp

gamma_VMn = mu_VMn*kb*T*log(eta)*(C_VMn_ox_g-C_VMn_ox_m*eta^ZVMn)/(1-eta^ZVMn);

kp = Vspinel*abs(gamma_VMn)    % Parabolic growth rate [nm²/hr]

%% Oxide thickness

delta = sqrt(2*kp*t)      %Oxide thickness after time t [nm]


%% Metal/oxide concentration

C_Mn_steel = x_steel_Mn* Csteel;                          % Native concentration of Mn atoms in steel [/nm^3]
DMn = 3600*10;                                            % Diffusion coefficient of Mn in 21-2N steel [nm^2/hr]
C_Mn_m_ox = C_Mn_steel - Cspinel*sqrt(pi/2*kp/DMn)