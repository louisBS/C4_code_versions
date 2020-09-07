%% C4_generate_MC
%%Generate the Monte Carlo transitions to run the C4 model on them

%%Author : Ryan Schulte, Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables

n_MC_trans = 1000;                                                   % Number of Monte Carlo transitions generated
d_total_MC = zeros(n_MC_trans, N+1);                                % Matrix containing the N+1 (time) total oxide thickness for the n_curves Monte Carlo transitions [cm]
std_gaussian = 0.15 * 1e-4;                                         % Standard deviation of the Gaussian distribution chosen to generate the transition thicknesses [cm]
mean_gaussian = 3.1 * 1e-4;                                         % Mean of the Gaussian distribution chosen to generate the transition thicknesses [cm]
num_transitions = 2;                                                % Number of transitions that are used for Monte Carlo study
transitions_MC = zeros(n_MC_trans, num_transitions);                % Matrix containing the n_curves Monte Carlo transitions thicknesses [cm] for the num_transisions transitions that are studied
time_trans = zeros(n_MC_trans, num_transitions);                    % Time step at which the transitions studied happen
d_diff_C4_exp_MC_total = zeros(n_MC_trans,1);                       % Normalized root mean square error (NRMSE) for each of the oxide thickness curves


%% Generate the transitions thicknesses following a Gaussian distribution

for c = 1:n_MC_trans      % Loop through the Monte Carlo transitions
    for trans = 1:num_transitions       % Loop through the number of transitions that are used for Monte Carlo study
        transitions_MC(c, trans) = normrnd(mean_gaussian, std_gaussian);    % The transition follows a Gaussian (or Normal) distribution
    end
end

%% Run the C4 model using the transition thickness generated as parameters

progress = waitbar(0, sprintf('run %d/%d', 0, n_MC_trans), 'Name', 'Monte Carlo generation ...');        % Initialize progress bar

for c = 1:n_MC_trans
    
    waitbar(c / n_MC_trans, progress, sprintf('run %d/%d', c, n_MC_trans));
    
    transition = ones(1, 4);
    for trans = 1:num_transitions
        transition(trans) = transitions_MC(c, trans);       % Generate the transition vector for each Monte Carlo run [cm] ; the "1" at the end is added artifically as a value too big to be reached
    end
    
    C4_parameter;                       % Definition of the constants and the parameters of the model /!\ The lines defining the transitions need to be commented to make this code run /!\
    C4_oxidation;                       % Main code creating the finite-element mesh and calling C4_ox_flux to solve the vacancy, electron and hydrogen fluxes in the oxide; saves the data in the output file
    
    d_total_MC(c, :) = d_total(:);      % Extract the total oxide thickness for every Monte Carlo run
    for trans = 1:num_transitions
        time_trans(c, trans) = N_transition(1, trans+1);   % Extract the time step at which the transitions happen
    end
    
    C4_exp_data;								% Load the experimental data associated with the alloy
    C4_least_square;                            % Least square criteria used to compare the C4 model with empirical models and experimental data
    d_diff_C4_exp_MC_total(c, 1) = d_diff_C4_exp;     % Extract the NRMSE computed by the C4_least_square code
end

close(progress);
clear progress;

save(strcat('output_file/', file_name, '_Gaussian'));           % Save the data in a separate file from the one used by the regular C4 code