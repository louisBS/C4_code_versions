%% C4_launch
%%Launch the different scripts of the C4 model

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

clearvars; clc; close;

%% Define the alloy, the model, the temperature and the duration of the oxidation

%%There are 3 alloys currently available. Each alloy has its own set of parameters (mainly migration energies) and experimental data.
%Zr05Nb is used for operating condition (360C, 240d), without any transition
%Zr10Nb is used for operating condition (360C, 240d), with several transitions
%Zry4 is used for high-temperature oxidation (950C/1100C/1200C, 1500s)

%%There are 4 different models, the difference between them is which effect is taken into account:
%C4-O takes into account the oxide dissolution into the metal -> necessary for high-temperature, not very useful for operating temperature (with exceptions).
%C4-H takes into account the transport of hydrogen through the oxide -> used for operating temperature if hydrogen content is needed, not very useful for high-temperature.
%C4-OH takes into account both effects of oxide dissolution into the metal and hydrogen pick-up
%C4 doesn't take into account neither the oxide dissolution, nor the hydrogen pick-up

%%There are different modes available:
%"Isotherm" is used for isothermal corrosion at operating or high temperature
%"Linear" is used for simple linear temperature transients ; the temperature string is the initial and final temperature separated by an hyphen
%"LOCA" is used for more complicated transients that try to reproduce a LOCA ; the different cases are defined in the variable temperature (Peak1200C ; CINOG5)
%"H optimization" is used for the optimization of the hydrogen migration energy at operating temperature
%"D HT optimization" is used for the optimization of the oxygen diffusion coefficient in the alpha and beta phase at high temperature
%"Monte Carlo" is used for the Monte Carlo transition study

%{
alloy = 'Zr05Nb';
model = 'C4-OH';
temperature = '360C';
exposure_time = '240d';
mode = 'Isotherm';
%}
%%{
alloy = 'Zry4';
model = 'C4-O';
temperature = '1200C';
exposure_time = '1500s';
mode = 'Isotherm';
%}
%{
alloy = 'Zry4';
model = 'C4-O';
temperature = 'LeistikowLOCA';
% exposure_time = '1800s';
mode = 'LOCA';
%}

% Prefix of the name of the file holding the input data and output data
if strncmp(temperature, 'CINOG', 5) || strncmp(temperature, 'Peak', 4) || strcmp(temperature, 'LeistikowLOCA')
    file_name = strcat(alloy, '_', model, '_', temperature);
else
    file_name = strcat(alloy, '_', model, '_', temperature, '_', exposure_time);
end

fprintf('                                                                                                        dddddddd                            \n');
fprintf('        CCCCCCCCCCCCC       444444444       MMMMMMMM               MMMMMMMM                             d::::::d                    lllllll \n');
fprintf('     CCC::::::::::::C      4::::::::4       M:::::::M             M:::::::M                             d::::::d                    l:::::l \n');
fprintf('   CC:::::::::::::::C     4:::::::::4       M::::::::M           M::::::::M                             d::::::d                    l:::::l \n');
fprintf('  C:::::CCCCCCCC::::C    4::::44::::4       M:::::::::M         M:::::::::M                             d:::::d                     l:::::l \n');
fprintf(' C:::::C       CCCCCC   4::::4 4::::4       M::::::::::M       M::::::::::M   ooooooooooo       ddddddddd:::::d     eeeeeeeeeeee     l::::l \n');
fprintf('C:::::C                4::::4  4::::4       M:::::::::::M     M:::::::::::M oo:::::::::::oo   dd::::::::::::::d   ee::::::::::::ee   l::::l \n');
fprintf('C:::::C               4::::4   4::::4       M:::::::M::::M   M::::M:::::::Mo:::::::::::::::o d::::::::::::::::d  e::::::eeeee:::::ee l::::l \n');
fprintf('C:::::C              4::::444444::::444     M::::::M M::::M M::::M M::::::Mo:::::ooooo:::::od:::::::ddddd:::::d e::::::e     e:::::e l::::l \n');
fprintf('C:::::C              4::::::::::::::::4     M::::::M  M::::M::::M  M::::::Mo::::o     o::::od::::::d    d:::::d e:::::::eeeee::::::e l::::l \n');
fprintf('C:::::C              4444444444:::::444     M::::::M   M:::::::M   M::::::Mo::::o     o::::od:::::d     d:::::d e:::::::::::::::::e  l::::l \n');
fprintf('C:::::C                        4::::4       M::::::M    M:::::M    M::::::Mo::::o     o::::od:::::d     d:::::d e::::::eeeeeeeeeee   l::::l \n');
fprintf(' C:::::C       CCCCCC          4::::4       M::::::M     MMMMM     M::::::Mo::::o     o::::od:::::d     d:::::d e:::::::e            l::::l \n');
fprintf('  C:::::CCCCCCCC::::C          4::::4       M::::::M               M::::::Mo:::::ooooo:::::od::::::ddddd::::::dde::::::::e          l::::::l\n');
fprintf('   CC:::::::::::::::C        44::::::44     M::::::M               M::::::Mo:::::::::::::::o d:::::::::::::::::d e::::::::eeeeeeee  l::::::l\n');
fprintf('     CCC::::::::::::C        4::::::::4     M::::::M               M::::::M oo:::::::::::oo   d:::::::::ddd::::d  ee:::::::::::::e  l::::::l\n');
fprintf('        CCCCCCCCCCCCC        4444444444     MMMMMMMM               MMMMMMMM   ooooooooooo      ddddddddd   ddddd    eeeeeeeeeeeeee  llllllll\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Launch the scripts
%%{
if strcmp(mode, 'Isotherm') || strcmp(mode, 'Linear') || strcmp(mode, 'LOCA')
    C4_parameter;	% Definition of the constants and the parameters of the model
    C4_oxidation;	% Main code creating the finite-element mesh and calling C4_ox_flux to solve the vacancy, electron and hydrogen fluxes in the oxide; saves the data in the output file
    
    load(strcat('output_file/', file_name));	% Make sure that the correct output file is opened for the analysis
    C4_exp_data;								% Load the experimental data associated with the alloy
    if strcmp(mode, 'Isotherm')
        %C4_least_square;							% Least square criteria used to compare the C4 model with empirical models and experimental data
        C4_MOOSE_comparison;
    end
    C4_plot;									% Generate all the plots
end
%}

%% Launch the scripts for the optimization of the hydrogen migration energy at operating temperature
%%{
if strcmp(mode, 'H optimization')
    C4_parameter;	% Definition of the constants and the parameters of the model
    C4_oxidation;	% Main code creating the finite-element mesh and calling C4_ox_flux to solve the vacancy, electron and hydrogen fluxes in the oxide; saves the data in the output file
    
    load(strcat('output_file/', file_name));	% Make sure that the correct output file is opened for the analysis
    C4_exp_data;								% Load the experimental data associated with the alloy
    
    C4_H_optimization;
    
    C4_plot;
end
%}


%% Launch the scripts for the optimization of the oxygen diffusion coefficient in the alpha and beta phase at high temperature
%%{
if strcmp(mode, 'D HT optimization')
    C4_parameter;	% Definition of the constants and the parameters of the model
    
    C4_D_HT_optimization;
    
    load(strcat('output_file/', file_name));	% Make sure that the correct output file is opened for the analysis
    C4_exp_data;								% Load the experimental data associated with the alloy
    C4_least_square;							% Least square criteria used to compare the C4 model with empirical models and experimental data
    C4_plot;									% Generate all the plots
end
%}


%% Launch the scripts for the Monte Carlo transition

%%{
if strcmp(mode, 'Monte Carlo')
    C4_parameter;   % Definition of the constants and the parameters of the model
%     C4_generate_MC; % Generate the Monte Carlo transitions to run the C4 model on
    
    load(strcat('output_file/', file_name,'_Gaussian'));
    C4_stats_MC;    % Statistical analysis of the Monte Carlo transitions
    C4_plot_MC;     % Generate all the plots linked to the Monte Carlo transitions
end
%}

toc