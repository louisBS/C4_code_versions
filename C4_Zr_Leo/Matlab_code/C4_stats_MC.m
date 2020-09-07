%% C4_stats_MC
%%Statistical analysis of the data generated by the Monte Carlo transitions code

%%Author : Ryan Schulte, Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Bin Sampling Method
%%Divide exposure time in time steps, total oxide thickness space into bins, and count the number of curves in each bin

dd_bin = .1 * 1e-4;                                                 % Width of the bins in the oxide thickness space [cm]
max_d_bins = 12 * 1e-4;                                             % Maximum oxide thickness used to define the bins [cm]
d_bins = 0:dd_bin:max_d_bins;                                       % oxide thickness vector using the oxide thickness bin [cm]
n_d_bins = length(d_bins);                                          % Number of oxide thickness bins
n_curve_per_bin = zeros(n_d_bins, N+1);                             % Matrix containing the number of curves that pass through each one of the n_d_bins oxide thickness bin at the N+1 time steps
bin_numbers = zeros(n_MC_trans, N+1);                               % Number of the oxide thickness bin in which each one of the n_MC_trans curves is at the N+1 time steps

% Assign every oxide thickness curve to an oxide thickness bin at each time step
for n = 1:N+1   % loop through exposure time
    bin_numbers(:,n) = discretize(d_total_MC(:, n), d_bins);     % discretize the oxide thickness space and assign the bin number corresponding to the oxide thickness at this time step
end

% Count the number of curves in each bin
for n = 1:N+1   % Loop through exposure time
    for c = 1:n_MC_trans    % Loop through oxide thickness curves
        for o = 1:n_d_bins    % Loop through the oxide thickness bins
            if o == bin_numbers(c, n)                                   % Checks if the curve is in the current bin
                n_curve_per_bin(o, n) = n_curve_per_bin(o, n) + 1 / n_MC_trans;      % Increment the count for the number of curves if the curve is in the current bin and normalize it to the total number of oxide thickness curves generated
            end
        end
    end
end


%% PDF (Probability Density Function) Fit Method
%%Fit a statistical distribution function to the distribution of oxide thicknesses (Number of oxide thickness curves that have this oxide thickness) at each time step

pdf_step = .01 * 1e-4;                                  % Size of the sampling bin for pdf [cm]
d_pdf = (0:pdf_step:max_d_bins)';                       % oxide thickness samplling for pdf [cm]
n_d_pdf = length(d_pdf);                                % Number of sampling bin for pdf
d_pdf_MC = zeros(n_d_pdf, N+1);                         % Oxide thickness distribution function
d_pdf_MC_transitioned = zeros(n_d_pdf, N+1);            % Oxide thickness distribution function of curves that already have transitioned
d_pdf_MC_non_transitioned = zeros(n_d_pdf, N+1);        % Oxide thickness distribution function of curves that have not yet transitioned

d_total_MC_transitioned = zeros(n_MC_trans, N+1);                   % During the transition time window, it will store the transitioned oxide thickness curves ; also used before the first transition
d_total_MC_non_transitioned = zeros(n_MC_trans, N+1);   % During the transition time window, it will store the non-transitioned oxide thickness curves

time_window_transition = zeros(num_transitions, 2);     % minimum and maximum time step at which the transition occurs for every transition
time_window = zeros(num_transitions, 2);            % time window considered for every transition ; the first line is for the pretransition

pdf_transitioned = cell(1, N+1);                        % Array containing the probability density function of the number of transitioned curves with each oxide thickness at every time step
pdf_non_transitioned = cell(1, N+1);                    % Array containing the probability density function of the number of non-transitioned curves with each oxide thickness at every time step

non_transitioned_proba = zeros(1, N+1);                 % Probability associated to the number of non-transitioned curves at every time step

% Pre-transition: all the oxide thickness curves are exactly the same
%Create a Dirac PDF at the oxide thickness value

for n = 1:min(time_trans(:,1))
    d_pre_transition = floor(d_total_MC(1, n) / pdf_step);
    if d_pre_transition == 0
        d_pdf_MC(1, n) = 1;
    else
        d_pdf_MC(d_pre_transition, n) = 1;
    end
end

% After the first transition

for trans_num = 1:num_transitions   % Loop through each transition that occurs
    
    % Minimum and maximum time step at which the transition occurs
    time_window_transition(trans_num, :) = [min(time_trans(:,trans_num))+1 max(time_trans(:,trans_num))];      
    
    % Time window of interest for each transition
    if trans_num == 1                   % For the first transition
        time_window(trans_num, 1) = time_window_transition(trans_num, 1);           % Current transition time window beginning     
    else
        time_window(trans_num, 1) = time_window_transition(trans_num - 1, 2) + 1;       % Use the previous transition time window end to determine the beginning of the current transition time window
    end
    if trans_num == num_transitions     % For the last transition
        time_window(trans_num, 2) = N + 1;
    else
        time_window(trans_num, 2) = time_window_transition(trans_num, 2);
    end
    
    
    % Remove non transitioned oxide thickness curves from d_total_MC_transitioned %%%%%???%%%%%
    for n = time_window(trans_num, 1):time_window(trans_num, 2)  % Loop through time within transition window
        for c = 1:n_MC_trans    % Loop through oxide thickness curves
            if n <= time_trans(c,trans_num)
                % The curve has not transitioned yet
                d_total_MC_non_transitioned(c, n) = d_total_MC(c, n);
                non_transitioned_proba(n) = non_transitioned_proba(n) + 1 / n_MC_trans;
            else
                d_total_MC_transitioned(c, n) = d_total_MC(c, n); 
            end
        end
    end
    
    
    % Fit Gaussian distribution to transitioned and non-transitioned curves
    for n = time_window(trans_num, 1):time_window(trans_num, 2)   % Loop through exposure time within the time window considered for this transition
        % Transitioned curves
        temp = nonzeros(d_total_MC_transitioned(:,n));
        if isempty(temp)
        elseif length(temp) == 1 %Solution for when there is only 1 transitioned curve
            d_pdf_MC_transitioned(floor(temp / pdf_step), n) = 1;
        else
            pdf_transitioned{n} = fitdist(temp, 'Normal');
        end
        % Sample the distribution
        if isempty(pdf_transitioned{n}) == 0
            if pdf_transitioned{n}.sigma < 1e-8
                d_pdf_MC_transitioned(floor(pdf_transitioned{n}.mu / pdf_step), n) = 1;
            else
                d_pdf_MC_transitioned(:, n) = pdf(pdf_transitioned{n}, d_pdf) * pdf_step; %Generates the distribution values for the probability density object 'pdArray' evaluated at the points within 'x_pdf' ; the normalization is over the integral.
                d_pdf_MC_transitioned(:, n) = d_pdf_MC_transitioned(:, n) ./ max(cumtrapz(d_pdf_MC_transitioned(:, n)));   % Renormalize to the cumulative integral = 1
            end
        end
        % non-transitioned curves
        temp = nonzeros(d_total_MC_non_transitioned(:,n));
        if isempty(temp)
        elseif length(temp) == 1 %Solution for when there is only 1 transitioned curve
            d_pdf_MC_non_transitioned(floor(temp / pdf_step), n) = 1;
        else
            pdf_non_transitioned{n} = fitdist(temp, 'Normal');
        end
        % Sample the distribution
        if isempty(pdf_non_transitioned{n}) == 0
            if pdf_non_transitioned{n}.sigma < 1e-8
                d_pdf_MC_non_transitioned(floor(pdf_non_transitioned{n}.mu / pdf_step), n) = 1;
            else
                d_pdf_MC_non_transitioned(:, n) = pdf(pdf_non_transitioned{n}, d_pdf) * pdf_step; %Generates the distrubtion values for the probability density object 'pdArray' evaluated at the points within 'x_pdf'
                d_pdf_MC_non_transitioned(:, n) = d_pdf_MC_non_transitioned(:, n) / max(cumtrapz(d_pdf_MC_non_transitioned(:, n)));   % Renormalize to the cumulative integral = 1
            end
        end
    end
end


%Re add the non-transitioned probability and scale each probability
d_pdf_MC = d_pdf_MC + (1 - non_transitioned_proba) .* d_pdf_MC_transitioned + non_transitioned_proba .* d_pdf_MC_non_transitioned;


%% Confidence Interval
%Strategy: Use cumulatie probabilities to find the region that contains a certain percentage of the curves 

confidence_bounds = zeros(2, N+1); % Lower bound in row 1 ; upper bound in row 2
d_cdf_MC = zeros(n_d_pdf, N+1); %Each column is the cumulative probability for an exposure time.

confidence_level = 0.90;

%Find the values that need to be checked for in the cumulative probabilities:
bound_threshold = [(1-confidence_level)/2, 1 - ((1-confidence_level)/2)];

for n = 1:N+1   % Loop through exposure times
    
    if n < min(time_trans(:,1)) % Pretransition region
        confidence_bounds(:, n) = [d_total_MC(1, n); d_total_MC(1, n)];
    else % After transitions start        
        d_cdf_MC(:, n) = cumtrapz(d_pdf_MC(:, n));
        
        %Finds the INDEX of the cumulative probabilities that are of the
        %appropriate size
        confidence_bounds(1, n) =  find(d_cdf_MC(:, n) > bound_threshold(1) * max(d_cdf_MC(:, n)), 1) * pdf_step;
        confidence_bounds(2, n) =  find(d_cdf_MC(:, n) > bound_threshold(2) * max(d_cdf_MC(:, n)), 1) * pdf_step;
    end
end


% Least square only on the time window of the transition
% Not really useful in fact ...

d_diff_C4_exp_MC_partial = zeros(n_MC_trans, num_transitions);      % Normalized root mean square error (NRMSE) for each of the oxide thickness curves

for trans_num = 1:num_transitions   % Loop through each transition that occurs
    day_wg_exp_partial = [];
    d_exp_partial = [];
    for k = 1:length(d_exp)
        if day_wg_exp(k) >= time_window(trans_num, 1) * dt / 86400 && day_wg_exp(k) < time_window(trans_num, 2) * dt / 86400
            day_wg_exp_partial = [day_wg_exp_partial day_wg_exp(k)];
            d_exp_partial = [d_exp_partial d_exp(k)];
        end
    end
    
    idx = zeros(1, length(d_exp_partial));
    for k = 1:length(d_exp_partial)
        t_curr = 86400 * day_wg_exp_partial(k);
        idx(k) = find(time == t_curr);
        for c = 1:n_MC_trans
            d_diff_C4_exp_MC_partial(c, trans_num)  = sqrt(d_diff_C4_exp_MC_partial(c, trans_num)^2 + ((1e4*d_total_MC(c, idx(k)) - d_exp_partial(k)) / d_exp_partial(k))^2 / length(d_exp_partial));
        end
    end
end

save(strcat('output_file/',file_name, '_Gaussian'));