%% C4_H_optimization
%%Fit hydrogen concentration by varying the hydrogen migration energy

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated: 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fit H concentration curve by changing Emh
%%{
progress = waitbar(0, sprintf('step %d/%d', 1, N + 1), 'Name', 'Optimizing ...');
i = 3;
optimization_factor = 0.1;
precision = 0.01;
max_intern_step = 20;

while i <= N + 1
    nb_intern_step = 0;
    while abs(Ch(i) - Ch_fit(i)) > precision && nb_intern_step < max_intern_step
        Emh(i-1:N + 1) = Emh(i-1:N + 1) + optimization_factor * (Ch(i) - Ch_fit(i));
        C4_oxidation;
        load(strcat('output_file/', file_name));
        C4_exp_data;
        nb_intern_step = nb_intern_step + 1;
    end
    
    if nb_intern_step == max_intern_step
        optimization_factor = optimization_factor / 2;
    else
        i = i + 1;
        optimization_factor = 0.1;
    end
    
    waitbar(i / (N + 1), progress, sprintf('step %d/%d', i, N + 1));
end

close(progress);

save(strcat('input_file/Emh_', alloy, '_', model, '_', temperature, '_', exposure_time), 'Emh');
%}

%% Fit H concentration curve by changing the hydrogen concentration at the oxide/water interface
%{
progress = waitbar(0, sprintf('step %d/%d', 1, N + 1), 'Name', 'Optimizing ...');
i = 3;
optimization_factor = 1e21;
precision = 0.1;
max_intern_step = 50;

while i <= N + 1
    nb_intern_step = 0;
    while abs(Ch(i) - Ch_fit(i)) > precision && nb_intern_step < max_intern_step
        Ch_ox_w(i-1:N + 1) = Ch_ox_w(i-1:N + 1) - optimization_factor * (Ch(i) - Ch_fit(i));
        C4_oxidation;
        load(strcat('output_file/', file_name));
        C4_exp_data;
        nb_intern_step = nb_intern_step + 1;
    end
    
    if nb_intern_step == max_intern_step
        optimization_factor = optimization_factor / 2;
    else
        i = i + 1;
        optimization_factor = 1e21;
    end
    
    waitbar(i / (N + 1), progress, sprintf('step %d/%d', i, N + 1));
end

Ch_ox_w_ppm = Ch_ox_w * Mh / (RhoZi * Na) * 1e6;        % Concentration of hydrogen in the oxide at the oxide/water interface [wt ppm]

close(progress);
toc

C4_plot;
%}

%% Fit both H concentration and oxide thickness curves by changing Eme and Emh in the same time
%{
for i = 4:N+1
    while abs(Ch(i)-Ch_fit(i)) > 0.1% || abs(wg_approx(i)-wg_fit(i)) > 1
        if abs(Ch(i)-Ch_fit(i)) > 0.1
            Emh(i:N+1) = Emh(i:N+1) + 0.001 * (Ch(i)-Ch_fit(i));
%             Eme(i:N+1) = Eme(i:N+1) - 0.001 * (Ch(i)-Ch_fit(i));

            
%             C4_ox_anisoth("C4");
%             C4_ox_anisoth("C4_H");
            C4_ox_anisoth("C4_OH");
            C4_exp_data;
            load(strcat('output_file/', file_name_C4_OH));
        end
%         if abs(wg_approx(i)-wg_fit(i)) > 1
%             Eme(i-1:N+1) = Eme(i-1:N+1) + 0.01 * (wg_approx(i)-wg_fit(i));
%             
% %             C4_ox_anisoth("C4");
%             C4_ox_anisoth("C4_H");
%             C4_plot;
%             load(strcat('output_file/', file_name_C4_H));
%         end
%         
% %         C4_ox_anisoth("C4");
%         C4_ox_anisoth("C4_H");
%         C4_plot;
%         load(strcat('output_file/', file_name_C4_H));
        
    end
end
%}