%% C4_plot
%%Generate the plots of the data

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set default graphic properties

set(0, 'DefaultFigureColor', 'White');
set(0, 'DefaultAxesFontSize', 16);      % 50 for presentation ; 30 else
set(0, 'DefaultAxesGridLineStyle', '-');
set(0, 'DefaultAxesXGrid', 'on');
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesXColor', 'Black');
set(0, 'DefaultAxesYColor', 'Black');
set(0, 'DefaultAxesZColor', 'Black');
set(0, 'DefaultAxesGridColor', 'Black');
set(0, 'DefaultAxeslineWidth', 1);
set(0, 'DefaultLineLineWidth', 2.5);      % 4 for presentation ; 3 else
set(0, 'DefaultLineMarkerSize', 10);    % 15 for presentation ; 10 else
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0.035 1 0.885]);
set(0, 'DefaultFigurePaperType', 'usletter');
set(0, 'DefaultFigurePaperOrientation', 'landscape');

% warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%%
file_name_C4 = strcat(alloy, '_C4_', temperature, '_', exposure_time);

if exposure_time(end) == 'd'
    time_plot = days;
    time_unit = 'days';
elseif exposure_time(end) == 's'
    time_plot = time;
    time_unit = 's';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Oxford TEM porosity counting
%{
fig_pore = figure();

distance = [300 900 1500];
pore_density = [2.57e-5 8.29e-5 1.5e-4];
std = [1.4e-5 1.67e-5 3.63e-5];

errorbar(distance, pore_density, std, '^k', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
title('Variation of porosity density in the oxide');
xlabel('Distance from the oxide-metal interface');
ylabel('Pore density (pore diameter / volume) [nm/nm^3]');
%}


%% Migration energies
%{
fig_Em = figure('Name', 'Migration energies');

hold on;
plot(time_plot, Emv, 'Color', [1 0.5 0], 'DisplayName', 'Vacancy migration energy');
plot(time_plot, Eme, 'b', 'DisplayName', 'Electron untrapping energy');
if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
    plot(time_plot, Emh, 'Color', [0 0.6 0], 'DisplayName', 'Hydrogen migration energy');
end

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, [1 2], '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylim([1 2]);
ylabel('Migration energy [eV]');
% title('Migration energies evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthEast');
%}


%% Diffusion coefficients in the oxide
%{
fig_D_ox = figure('Name', 'Diffusion coefficient in the oxide');

semilogy(time_plot, D_v * ones(1,N+1) * 1e-4, 'Color', [1 0.5 0], 'DisplayName', 'Vacancy');
hold on;
semilogy(time_plot, D_e * ones(1,N+1) * 1e-4, 'b', 'DisplayName', 'Electron');
if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
    semilogy(time_plot, D_h * 1e-4, 'Color', [0 0.6 0], 'DisplayName', 'Hydrogen');
    semilogy(time_plot, 6e-19 * ones(1, N+1), ':', 'Color', [0 0.6 0], 'DisplayName', 'Literature Zry-4 (Hatano)');
    semilogy(time_plot, 1.8e-19 * ones(1, N+1), '--', 'Color', [0 0.6 0], 'DisplayName', 'Literature Zr-2.5Nb (McIntyre)');
    semilogy(time_plot, 1.13e-17 * ones(1, N+1), '-.', 'Color', [0 0.6 0], 'DisplayName', 'Literature Zr-2.5Nb (Khatamian)');
    semilogy(time_plot, 1.42e-16 * ones(1, N+1), '-.', 'Color', [0 0.6 0], 'DisplayName', 'Literature Zr-2.5Nb (Khatamian)');
end

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, [1 2], '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylim([1e-19 1e-15]);
ylabel('Diffusion coefficient [m^2/s]');
% title('Diffusion coefficient in the oxide evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthWest');
%}


%% Temperature evolution
%{

if strcmp(mode, 'Linear') || strcmp(mode, 'LOCA')
    fig_T = figure('Name', 'Temperature');
    
    hold on;
    plot(time_plot, T - Tk, '-k', 'DisplayName', 'Temperature');
    
    xlim([0 time_plot(end)]);
    xlabel(strcat('Exposure time [', time_unit, ']'));
%     ylim([1 2]);
    ylabel('Temperature [C]');
%     title('Temperature evolution with exposure time');
%     legend('-DynamicLegend', 'Location', 'NorthEast');
end
%}


%% Particle flux
%{
fig_flux = figure('Name', 'Particle fluxes');

subplot(2,1,1);

semilogy(time_plot(2:end), Jv(2:end), 'Color', [1 0.5 0], 'DisplayName', 'Vacancy flux');
hold on;
semilogy(time_plot(2:end), Je(2:end), 'b', 'DisplayName', 'Electron flux');
if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
    semilogy(time_plot(2:end), Jh(2:end), 'Color', [0 0.6 0], 'DisplayName', 'Hydrogen flux');
end

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('Particle flux [cm^{-2}s^{-1}]');
title('Particle fluxes evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthEast');

subplot(2,1,2);
hold on;
plot(time_plot(2:end), 2*Jv(2:end) - Je(2:end) - Jh(2:end), 'k', 'DisplayName', 'Net current in the oxide');

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)])
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('Particle flux [cm^{-2}s^{-1}]');
title('Coupled current equation evolution with exposure time');
legend('-DynamicLegend', 'Location', 'SouthEast');
%}


%% Weight gain
%{
fig_wg = figure('Name', 'Weight gain');

hold on;
if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
    % High-temperature experimental data
    if exist('wg_exp_NRC', 'var') ==1
        plot(time_exp_NRC, wg_exp_NRC, '^', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'DisplayName', 'NRC experiment');
    end
    plot(time_exp_CP, wg_exp_CP, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Cathcart-Pawel experiment');
    if exist('wg_exp_UW', 'var') == 1
        plot(time_exp_UW, wg_exp_UW, '^r', 'MarkerFaceColor', 'r', 'DisplayName', 'UW experiment');
    end
    if exist('wg_exp_MIT', 'var') == 1
%         plot(time_exp_MIT, wg_exp_MIT, '^b', 'MarkerFaceColor', 'b', 'DisplayName', 'MIT experiment');
    end
    if exist('wg_exp_WPI', 'var') == 1
        plot(time_exp_WPI, wg_exp_WPI, '^g', 'MarkerFaceColor', 'g', 'DisplayName', 'Biederman (WPI) experiment');
    end
%     if exist('wg_exp_Brachet_Zry4', 'var') == 1
%         plot(time_wg_exp_Brachet_Zry4, wg_exp_Brachet_Zry4, '^c', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
%         plot(time_wg_exp_Brachet_M5, wg_exp_Brachet_M5, 'vc', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
%     end
%     if exist('wg_exp_ORNL', 'var') == 1
%         plot(time_exp_ORNL, wg_exp_ORNL, '^', 'Color', [1 0.1 0.6], 'MarkerFaceColor', [1 0.1 0.6], 'DisplayName', 'ORNL experiment');
%     end
%     plot(time_exp_Lemmon(1), wg_exp_Lemmon(1), '^m', 'MarkerFaceColor', 'm', 'DisplayName', 'Lemmon experiment');
    
    % High-temperature models
    plot(time_CP, wg_CP, '-k', 'DisplayName', 'Cathcart-Pawel model');
    CP_extended = plot(time_CP_extended, wg_CP_extended, '--k', 'DisplayName', 'Cathcart-Pawel model');
    set(get(get(CP_extended,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    plot(time_BJ, wg_BJ, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model');
    plot(time_Leistikow, wg_Leistikow, 'Color', [0.5 0 0.5], 'DisplayName', 'Leistikow model');
    plot(time_WPI, wg_WPI, 'g', 'DisplayName', 'Biederman model');
    plot(time_Kawasaki, wg_Kawasaki, '-', 'Color', [0 0.6 0], 'DisplayName', 'Kawasaki model');
    plot(time_Urbanic, wg_Urbanic, '-', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Urbanic model');
%     plot(time_Lemmon, wg_Lemmon, '-m', 'DisplayName', 'Lemmon model');
%     plot(time, wg_Parsons, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Parsons model');
%     plot(time, wg_Hobson, '-c', 'DisplayName', 'Hobson model');
%     plot(time, wg_Klepfer, '-b', 'DisplayName', 'Klepfer model');
    
    % BISON
    %{
    plot(t_bison_Cathcart, wg_bison_Cathcart * 1e2, 'm', 'DisplayName', 'Bison Cathcart-Pawel model');
    plot(t_bison_Leistikow, wg_bison_Leistikow * 1e2, '-b', 'DisplayName', 'Leistikow model');
    plot(t_bison_isotherm, wg_bison_isotherm * 1e2, 'c', 'DisplayName', 'Bison C4 isotherm model');
    plot(t_bison_C4_O, wg_bison_C4_O * 1e2, '--b', 'DisplayName', 'Bison C4-O model');
    plot(t_bison_C4_O_fine, wg_bison_C4_O_fine * 1e2, '-r', 'DisplayName', 'C4 model');
    %}
    
    elseif strcmp(temperature, '360C')
        % Operating temperature experimental data
        plot(day_wg_exp, wg_exp * 1e-2, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Experiment');
        plot(time_plot, wg_approx, ':r', 'DisplayName', 'C4 model (approximate)');
    end
    
elseif strcmp(mode, 'Linear')
    % Linear temperature transient experimental data
    if strcmp(alloy, 'Zry4') && strcmp(temperature, '700-1200C')
        plot(time_plot, wg_CINOG, 'b', 'DisplayName', 'CINOG Experiment');
        plot(time_plot, wg_CP, '-k', 'DisplayName', 'Cathcart-Pawel model');
        plot(time_plot, wg_BJ, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model');
    end
    
elseif strcmp(mode, 'LOCA')
    if strcmp(alloy, 'Zry4') && strncmp(temperature, 'CINOG', 5)        % CINOG data
        plot(time_plot(N+1), wg_exp_CINOG, '^b', 'MarkerFaceColor', 'b', 'DisplayName', 'Experiment');
        plot(time_plot(N+1), wg_CINOG, 'vb', 'MarkerFaceColor', 'b', 'DisplayName', 'CINOG');
    elseif strcmp(temperature, '1100C') || strcmp(temperature, '1200C') % ANL data
        plot(time_plot(N+1), wg_exp_ANL, '^b', 'MarkerFaceColor', 'b', 'DisplayName', 'ANL Experiment');
        plot(time_plot(N+1), wg_CP_ANL, 'vb', 'MarkerFaceColor', 'b', 'DisplayName', 'ANL CP prediction');
    elseif strncmp(temperature, 'Peak', 4)                              % Leistikow data
        plot(time_plot(N+1) * ones(1, length(wg_exp_Leistikow)), wg_exp_Leistikow, '^b', 'MarkerFaceColor', 'b', 'DisplayName', 'Experiment');
    end
    if strcmp(alloy, 'Zry4') && strcmp(temperature, 'LeistikowLOCA')
        plot(time_exp_Leistikow_LOCA, wg_exp_Leistikow_LOCA, '^b', 'MarkerFaceColor', 'b', 'LineWidth', 2, 'MarkerSize', 15, 'DisplayName', 'Leistikow LOCA experiment');
        plot(time_Leistikow_LOCA, wg_Leistikow_LOCA, 'b', 'DisplayName', 'Leistikow SIMTRAN Model');
    end
    plot(time_plot, wg_CP, '-k', 'DisplayName', 'Cathcart-Pawel model');
    plot(time_plot, wg_BJ, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model');
end


% C4 model
plot(time_plot, wg, 'r', 'DisplayName', 'C4 model');
% % % plot(time_plot, wg_int, '--r', 'DisplayName', 'C4 model integral');
% plot(time_plot, wg_fit, '--r', 'DisplayName', 'C4 model power fit');

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
% if exist('wg_BJ', 'var') ==1
%     ylim([0 wg_BJ(length(time))]);
% else
%     ylim([0 inf]);
% end
ylabel('Weight gain [mg/cm^2]');
% title('Weight gain evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthWest');

% print(fig_wg, 'Zry-4_wg_1200C_1500s', '-depsc');
%}


%% Weight gain per sqrt(time)
%{
fig_wg_sqrt_t = figure('Name', 'Weight gain vs sqrt(time)');

hold on;
if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
    % High-temperature experimental data
    if exist('wg_exp_NRC', 'var') ==1
        plot(sqrt(time_exp_NRC), wg_exp_NRC, '^', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'DisplayName', 'NRC experiment');
    end
    plot(sqrt(time_exp_CP), wg_exp_CP, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Cathcart-Pawel experiment');
    if exist('wg_exp_UW', 'var') == 1
        plot(sqrt(time_exp_UW), wg_exp_UW, '^r', 'MarkerFaceColor', 'r', 'DisplayName', 'UW experiment');
    end
    if exist('wg_exp_WPI', 'var') == 1
        plot(sqrt(time_exp_WPI), wg_exp_WPI, '^g', 'MarkerFaceColor', 'g', 'DisplayName', 'Biederman (WPI) experiment');
    end
    if exist('wg_exp_Brachet_Zry4', 'var') == 1
        plot(sqrt(time_wg_exp_Brachet_Zry4), wg_exp_Brachet_Zry4, '^c', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
        plot(sqrt(time_wg_exp_Brachet_M5), wg_exp_Brachet_M5, 'vc', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
    end
    if exist('wg_exp_ORNL', 'var') == 1
        plot(sqrt(time_exp_ORNL), wg_exp_ORNL, '^', 'Color', [1 0.1 0.6], 'MarkerFaceColor', [1 0.1 0.6], 'DisplayName', 'ORNL experiment');
    end
%     plot(sqrt(time_exp_Lemmon(1)), wg_exp_Lemmon(1), '^m', 'MarkerFaceColor', 'm', 'DisplayName', 'Lemmon experiment');
    
    % High-temperature models
    plot(sqrt(time_CP), wg_CP, '-k', 'DisplayName', 'Cathcart-Pawel model');
    CP_extended = plot(sqrt(time_CP_extended), wg_CP_extended, '--k', 'DisplayName', 'Cathcart-Pawel model');
    set(get(get(CP_extended,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    plot(sqrt(time_BJ), wg_BJ, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model');
    plot(sqrt(time_Leistikow), wg_Leistikow, 'Color', [0.5 0 0.5], 'DisplayName', 'Leistikow model');
    plot(sqrt(time_WPI), wg_WPI, 'g', 'DisplayName', 'Biederman model');
    plot(sqrt(time_Kawasaki), wg_Kawasaki, '-', 'Color', [0 0.6 0], 'DisplayName', 'Kawasaki model');
    plot(sqrt(time_Urbanic), wg_Urbanic, '-', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Urbanic model');
%     plot(sqrt(time_Lemmon), wg_Lemmon, '-m', 'DisplayName', 'Lemmon model');
%     plot(sqrt(time_plot), wg_Parsons, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Parsons model');
%     plot(sqrt(time_plot), wg_Hobson, '-c', 'DisplayName', 'Hobson model');
%     plot(sqrt(time_plot), wg_Klepfer, '-b', 'DisplayName', 'Klepfer model');
    
    % BISON
    %{
    plot(sqrt(t_bison_Cathcart), wg_bison_Cathcart * 1e2, 'm', 'DisplayName', 'Bison Cathcart-Pawel model');
    plot(sqrt(t_bison_Leistikow), wg_bison_Leistikow * 1e2, '-b', 'DisplayName', 'Leistikow model');
    plot(sqrt(t_bison_isotherm), wg_bison_isotherm * 1e2, 'c', 'DisplayName', 'Bison C4 isotherm model');
    plot(sqrt(t_bison_C4_O), wg_bison_C4_O * 1e2, '--b', 'DisplayName', 'Bison C4-O model');
    plot(sqrt(t_bison_C4_O_fine), wg_bison_C4_O_fine * 1e2, '-r', 'DisplayName', 'C4 model');
    %}
    elseif strcmp(temperature, '360C')
        % Operating temperature experimental data
        plot(sqrt(day_wg_exp), wg_exp * 1e-2, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Experiment');
        plot(sqrt(time_plot), wg_approx, ':r', 'DisplayName', 'C4 model (approximate)');
    end
elseif strcmp(mode, 'Linear')
    % Temperature transient experimental data
    if strcmp(alloy, 'Zry4') && strcmp(temperature, '700-1200C')
        plot(sqrt(time_plot), wg_CINOG, 'b', 'DisplayName', 'CINOG Experiment');
    end
end


% C4 model
plot(sqrt(time_plot), wg, 'r', 'DisplayName', 'C4 model');
plot(sqrt(time_plot), wg_int, '--r', 'DisplayName', 'C4 model integral');
% plot(sqrt(time_plot), wg_fit, '--r', 'DisplayName', 'C4 model power fit');

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 sqrt(time_plot(end))]);
xlabel(strcat('Square root of exposure time [', time_unit, '^{1/2}]'));
ylabel('Weight gain [mg/cm^2]');
title('Weight gain evolution with square root of exposure time');
legend('-DynamicLegend', 'Location', 'SouthEast');
%}


%% Oxide thickness
%{
fig_d = figure('Name', 'Oxide thickness');

hold on;
if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
        % High-temperature experimental data
        if exist('d_exp_NRC', 'var') == 1
            plot(time_exp_NRC, d_exp_NRC, '^', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'DisplayName', 'NRC experiment');
        end
        plot(time_exp_CP, d_exp_CP, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Cathcart-Pawel experiment');
        if exist('d_exp_UW', 'var') == 1
            plot(time_exp_UW, d_exp_UW, '^r', 'MarkerFaceColor', 'r', 'DisplayName', 'UW experiment');
        end
        %         if exist('d_exp_MIT', 'var') == 1
        %             plot(time_exp_MIT, d_exp_MIT, '^b', 'MarkerFaceColor', 'b', 'DisplayName', 'MIT experiment');
        %         end
        %         if exist('d_exp_Brachet_Zry4', 'var') == 1
        %             plot(time_d_exp_Brachet_Zry4, d_exp_Brachet_Zry4, '^c', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
        %             plot(time_d_exp_Brachet_M5, d_exp_Brachet_M5, 'vc', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
        %         end
        %         if exist('d_exp_ORNL', 'var') == 1
        %             plot(time_exp_ORNL, d_exp_ORNL, '^', 'Color', [1 0.1 0.6], 'MarkerFaceColor', [1 0.1 0.6], 'DisplayName', 'ORNL experiment');
        %         end
        %         plot(time_exp_Lemmon(1), d_exp_Lemmon(1), '^m', 'MarkerFaceColor', 'm', 'DisplayName', 'Lemmon experiment');
        %         plot(time_exp_Sawarn, 0.5*d_exp_Sawarn, '^', 'Color', [1 0.1 0.58], 'MarkerFaceColor', [1 0.1 0.58], 'DisplayName', 'Sawarn experiment');
        
        % High-temperature models
        plot(time_CP, d_CP * 1e4, 'k', 'DisplayName', 'Cathcart-Pawel model');
        CP_extended = plot(time_CP_extended, d_CP_extended * 1e4, '--k', 'DisplayName', 'Cathcart-Pawel model');
        set(get(get(CP_extended,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        plot(time_BJ, d_BJ * 1e4, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model');
        plot(time_Leistikow, d_Leistikow * 1e4, 'Color', [0.5 0 0.5], 'DisplayName', 'Leistikow model');
        plot(time_WPI, d_WPI * 1e4, 'g', 'DisplayName', 'Biederman model');
        plot(time_Kawasaki, d_Kawasaki * 1e4, 'Color', [0 0.6 0], 'DisplayName', 'Kawasaki model');
        plot(time_Urbanic, d_Urbanic * 1e4, '-', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Urbanic model');
        
        % BISON
        %{
        plot(t_bison_Cathcart, d_bison_Cathcart * 1e6, '--g', 'DisplayName', 'Cathcart-Pawel model (BISON)');
%         plot(t_bison_Leistikow, d_bison_Leistikow * 1e6, '-b', 'DisplayName', 'Leistokow model');
%         plot(t_bison_isotherm, d_bison_isotherm * 1e6, 'c', 'DisplayName', 'Bison C4 isotherm model');
        plot(t_bison_C4_O, d_bison_C4_O * 1e6, '--b', 'DisplayName', 'C4 model (BISON)');
%         plot(t_bison_C4_O_fine, d_bison_C4_O_fine * 1e6, '-r', 'DisplayName', 'C4 model');
        %}
        
    elseif strcmp(temperature, '360C')
        % Operating temperature experimental data
        plot(day_wg_exp, d_exp, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Experiment');
%         plot(time_plot, wg / 0.1477, ':r', 'DisplayName', 'C4 model');
    end
    
    
elseif strcmp(mode, 'Linear') || strcmp(mode, 'LOCA')
    if strcmp(alloy, 'Zry4') && strcmp(temperature, 'LeistikowLOCA')
        errorbar(202, d_exp_Leistikow_LOCA, 2, '^b', 'MarkerFaceColor', 'b', 'LineWidth', 2, 'MarkerSize', 15, 'DisplayName', 'Leistikow LOCA experiment');
        plot(time_Leistikow_LOCA, d_Leistikow_LOCA, 'b', 'MarkerFaceColor', 'r', 'DisplayName', 'Leistikow SIMTRAN Model');
    end
    % Empirical models
    plot(time_plot, d_CP * 1e4, '-k', 'DisplayName', 'Cathcart-Pawel model');
    plot(time_plot, d_BJ * 1e4, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model');
end


% C4 model
plot(time_plot, d_total * 1e4, 'r', 'DisplayName', 'C4 model');
% plot(time_plot, d_fit * 1e4, '--r', 'DisplayName', 'C4 model power fit');

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end


xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('Oxide thickness [\mum]');
% title('Oxide thickness evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthWest');

% print(fig_d, 'Zry-4_d_1200C_1500s', '-depsc');
%}


%% Alpha-phase layer thickness
%{
fig_alpha = figure('Name', 'Alpha-phase layer thickness');

hold on;
if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
        % High-temperature experimental data
        if exist('alpha_exp_NRC', 'var') == 1
            plot(time_exp_NRC, alpha_exp_NRC, '^', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'DisplayName', 'NRC experiment');
        end
        plot(time_exp_CP, alpha_exp_CP, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Cathcart-Pawel experiment');
        if exist('alpha_exp_UW', 'var') == 1
            plot(time_exp_UW, alpha_exp_UW, '^r', 'MarkerFaceColor', 'r', 'DisplayName', 'UW experiment');
        end
%         if exist('alpha_exp_MIT', 'var') == 1
%             plot(time_exp_MIT, alpha_exp_MIT, '^b', 'MarkerFaceColor', 'b', 'DisplayName', 'MIT experiment');
%         end
%         if exist('alpha_exp_Brachet_Zry4', 'var') == 1
%             plot(time_alpha_exp_Brachet_Zry4, alpha_exp_Brachet_Zry4, '^c', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
%             plot(time_alpha_exp_Brachet_M5, alpha_exp_Brachet_M5, 'vc', 'MarkerFaceColor', 'c', 'DisplayName', 'Brachet experiment');
%         end
%         plot(time_exp_Sawarn, 0.5*alpha_exp_Sawarn, '^', 'Color', [1 0.1 0.58], 'MarkerFaceColor', [1 0.1 0.58], 'DisplayName', 'Sawarn experiment');
        
        % High-temperature models
        plot(time_CP, alpha_CP * 1e4, 'k', 'DisplayName', 'Cathcart-Pawel model');
        CP_extended = plot(time_CP_extended, alpha_CP_extended * 1e4, '--k', 'DisplayName', 'Cathcart-Pawel model');
        set(get(get(CP_extended,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        plot(time_WPI, alpha_WPI * 1e4, 'g', 'DisplayName', 'Biederman model');
        plot(time_Kawasaki, alpha_Kawasaki * 1e4, 'Color', [0 0.6 0], 'DisplayName', 'Kawasaki model');
        plot(time_Urbanic, alpha_Urbanic * 1e4, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Urbanic model');
    end

elseif strcmp(mode, 'Linear') || strcmp(mode, 'LOCA')
    if strcmp(alloy, 'Zry4') && strcmp(temperature, 'LeistikowLOCA')
        errorbar(202, alpha_exp_Leistikow_LOCA, 2, '^b', 'MarkerFaceColor', 'b', 'LineWidth', 2, 'MarkerSize', 15, 'DisplayName', 'Leistikow LOCA experiment');
        plot(time_Leistikow_LOCA, alpha_Leistikow_LOCA, 'b', 'MarkerFaceColor', 'r', 'DisplayName', 'Leistikow SIMTRAN Model');
    end
    % Empirical models
    plot(time_plot, d_CP * 1e4, '-k', 'DisplayName', 'Cathcart-Pawel model');
    plot(time_plot, d_BJ * 1e4, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model');
end

% C4 model
plot(time_plot, alpha * 1e4, 'r', 'DisplayName', 'C4 model');
% plot(time_plot, alpha_fit * 1e4, '--r', 'DisplayName', 'C4 model power fit');

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('\alpha-phase layer thickness [\mum]');
% title('\alpha-phase layer thickness evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthWest');

% print(fig_alpha, 'Zry-4_alpha_1200C_1500s', '-depsc');
%}


%% Ductile beta layer thickness
%{
fig_ductile_beta = figure('Name', 'Ductile beta layer thickness');

hold on;
plot(time_plot, 100 * duct_b / e_Zr, 'r', 'DisplayName', 'C4 model');

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('ductile \beta-phase thickness ratio [%]');
title('ductile \beta-phase thickness ratio evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthEast');
%}


%% Electric field
%{
fig_E = figure('Name', 'Electric field');

hold on;
plot(time_plot(2:end), E(2:end), 'r', 'DisplayName', strcat(model, ' model'));

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('Electric field [V/cm]');
title('Electric field evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthWest');
%}


%% Ratio of fluxes with 2 driving forces over Fickian flux
%{
fig_flux_ratio = figure('Name', 'Ratio of electrochemical flux over Fickian flux');

D_v = mu_v * kb * T_avg / Zv;
D_e = mu_e * kb * T_avg / Ze;

subplot(2,1,1);
semilogy(time_plot(2:end), Jv(2:end), 'Color', 'r', 'DisplayName', 'Total vacancy flux');
hold on;
semilogy(time_plot(2:end), D_v * (Cv_ox_a - Cv_ox_w) ./ d_protect(2:end), 'Color', 'b', 'DisplayName', 'Fickian vacancy flux');

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('Flux [cm^{-2}s^{-1}]');
title('Total and fickian flux');
legend('-DynamicLegend', 'Location', 'NorthWest');


subplot(2,1,2);
plot(time_plot(2:end), Jv(2:end) ./ (D_v * (Cv_ox_a - Cv_ox_w) ./ d_protect(2:end)), 'Color', 'r', 'DisplayName', 'Vacancy flux ratio');
hold on;
% semilogy(time_plot(2:end), Je(2:end) ./ (D_e * (Ce_ox_a - Ce_ox_w) ./ d_protect(2:end)), 'b', 'DisplayName', 'Electron flux ratio');
% if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
%     D_h = mu_h * kb * T_avg / Zh;
%     semilogy(time_plot(2:end), D_h * (Ch_ox_w - Ch_ox_a) ./ d_protect(2:end), 'Color', [0 0.6 0], 'DisplayName', 'Hydrogen flux ratio');
% end

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
title('Ratio of the total flux over the Fickian flux');
legend('-DynamicLegend', 'Location', 'NorthWest');
xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('');
%}


%% Hydrogen content
%{
if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
    fig_Ch = figure('Name', 'Hydrogen content');
    
    hold on;
    plot(day_Ch_exp, Ch_exp, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Experiment');
    plot(time_plot, Ch_fit, '--k', 'DisplayName', 'Experimental fit')
    plot(time_plot, Ch, 'r', 'DisplayName', 'C4 model');%strcat(model, ' model'));
    
    if length(transition) > 1
        transition_h = zeros(1,length(N_transition));
        for k = 2:length(N_transition)
            transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
            set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
    end
    
    xlim([0 time_plot(end)]);
    xlabel(strcat('Exposure time [', time_unit, ']'));
    ylabel('Hydrogen content [wt ppm]');
%     title('Hydrogen content evolution with exposure time');
    legend('-DynamicLegend', 'Location', 'NorthWest');
end
%}


%% Instantaneous Hydrogen pickup fraction
%{
if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
    fig_fh_inst = figure('Name', 'Instantaneous hydrogen pickup fraction');
    
    plot(time_plot, fh_inst * 1e2, 'r', 'DisplayName', strcat(model, ' model'));
    hold on;
    
    if length(transition) > 1
        transition_h = zeros(1,length(N_transition));
        for k = 2:length(N_transition)
            transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
            set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
    end
    
    xlim([0 time_plot(end)]);
    xlabel(strcat('Exposure time [', time_unit, ']'));
    ylabel('Hydrogen pickup fraction [%]');
    title('Instantaneous Hydrogen pickup fraction evolution with exposure time');
    legend('-DynamicLegend', 'Location', 'NorthWest');
end
%}


%% Total Hydrogen pickup fraction
%{
if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
    fig_fh_tot = figure('Name', 'Total hydrogen pickup fraction');
    
    errorbar(day_fh_exp, fh_exp*1e2, fh_exp_err*1e2, '^k', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'DisplayName', 'Experiment');
    hold on;
    plot(time_plot(2:end), fh_tot(2:end)*1e2, 'r', 'DisplayName', 'C4 model');
%     plot(time_plot, fh_fit*1e2, 'k', 'DisplayName', 'Experimental fit');
    
    if length(transition) > 1
        transition_h = zeros(1,length(N_transition));
        for k = 2:length(N_transition)
            transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
            set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
    end
    
    xlim([0 time_plot(end)]);
    xlabel(strcat('Exposure time [', time_unit, ']'));
    ylabel('Hydrogen pickup fraction [%]');
%     title('Total Hydrogen pickup fraction evolution with exposure time');
    legend('-DynamicLegend', 'Location', 'NorthWest');
end
%}


%% Evolution of hydrogen concentration at the oxide/water interface in case it is used as a fitting parameter
%{
fig_Ch_ox_w = figure('Name', 'Fit using Hydrogen concentration at the oxide/water interface');


semilogy(time_plot, Ch_ox_w_ppm, 'r', 'DisplayName', 'Vacancy migration energy');
hold on;

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, [5e2 2e4], '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
% ylim([1 2]);
ylabel('Hydrogen concentration [wt ppm]');
% title('Migration energies evolution with exposure time');
% legend('-DynamicLegend', 'Location', 'NorthEast');
%}


%% Difference between weight gain computed from the oxygen flux and from integration of concentration
%{
fig_wg_error = figure('Name', 'Weight gain difference');

wg_error = zeros(1, N+1);
for n = 2:N+1
    wg_error(1,n) = (wg(n)-wg_int(n)) - (wg(n-1)-wg_int(n-1));
end

hold on;
yyaxis left;
plot(time_plot, wg-wg_int,'r', 'DisplayName', 'Weight gain difference');


ylabel('Weight gain difference [mg/cm^2]');
set(gca, 'YColor', 'r');
set(gca, 'GridColor', 'k');

yyaxis right;
plot(time_plot, wg_error, 'b', 'DisplayName', 'Weight gain difference per step');

if length(transition) > 1
    transition_h = zeros(1,length(N_transition));
    for k = 2:length(N_transition)
        transition_h(k) = plot([N_transition(k)-1 N_transition(k)-1] * dt / 86400, ylim, '-.k', 'Color', [0.5 0.5 0.5], 'DisplayName', 'transition');
        set(get(get(transition_h(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    set(get(get(transition_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
end

xlim([0 time_plot(end)]);
xlabel(strcat('Exposure time [', time_unit, ']'));
% ylim([-10*wg_error(end) 10*wg_error(end)]);
ylabel('Weight gain difference [mg/cm^2]');
set(gca, 'YColor', 'b');

title('Weight gain difference evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthWest');
%}


%% 2-D Oxygen concentration profile
%%{
if strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')
    fig_Co = figure('Name', 'Oxygen concentration profile');
    
    n_time_plot = 3;
    plot_xo_interface = zeros(1,n_time_plot);
    plot_xo_ox = zeros(1,n_time_plot);
    plot_x = zeros(1,n_time_plot);
    list_colors = {'r', 'b', '0 0.6 0', [0.5 0.5 0.5], [1 0.5 0]};
%     t_plots = floor(linspace(1,N+1,n_time_plot));
    t_plots = [31 101 1501];
    
    %res_plot(1:K+2) = 1 ./ (1 + RhoZr * Na ./ (MZr .* res(1:K+2)));          % Convert oxygen concentration into oxygen fraction
    %res_plot(1:K+2) = res(1:K+2)/Czr;                             % convert O concentration to Co/Czr (strong discontinuity matlab variable)

    hold on;
    for k = 1:n_time_plot
        plot(1e4*x(:, t_plots(k)), 100*xo(:, t_plots(k))./(1-xo(:, t_plots(k))), 'Color', list_colors{k}, 'DisplayName', strcat(sprintf('t = %.1f ', time_plot(t_plots(k))), time_unit));
%         plot(1e4*x(:, t_plots(k)-1), 100*res_plot(1:K+2, t_plots(k)), ':', 'Color', list_colors{k}, 'DisplayName', strcat(sprintf('t = %.1f ', time_plot(t_plots(k))), time_unit));
%         plot_x(k) = plot(1e4*x(:,t_plots(k)), zeros(1,K+2), '+', 'Color', list_colors{k}, 'MarkerSize', 10, 'DisplayName', strcat(sprintf('t = %.1f ', time_plot(t_plots(k))), time_unit));
        plot_xo_interface(k) = plot([xi_ox(t_plots(k))*1e4 xi_ox(t_plots(k))*1e4], 100*[xo_a_ox xo_ox(1, t_plots(k))], 'Color', list_colors{k});
        plot_xo_ox(k) = plot(1e4*(xi_ox(t_plots(k))+x_ox(:, t_plots(k)))', 100*xo_ox(:, t_plots(k))', 'Color', list_colors{k});
        set(get(get(plot_xo_interface(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(plot_xo_ox(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%         set(get(get(plot_x(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    
    xlim([300 1e4*(e_Zr+d_total(end))]);
    xlabel('Cladding thickness [\mum]');
    ylim([0 100]);
    ylabel('Oxygen fraction [%]');
    title('Oxygen fraction profiles for different times');
    legend('-DynamicLegend', 'Location', 'NorthWest');
end
%}


%% 3-D Oxygen concentration profile
%{
if strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')
    fig_Co = figure('Name', '3-D Oxygen concentration profile');
    
    time_surf = zeros(1, N/10+1);
    x_surf = zeros(K+2+K_ox+1, N/10+1);
    xo_surf = zeros(K+2+K_ox+1, N/10+1);
    x_ox_surf = zeros(K_ox+1, N/10+1);
    xo_ox_surf = zeros(K_ox+1, N/10+1);
    
    for n = 1:N+1
        if mod(n,10) == 1
            time_surf(1,(n-1)/10+1) = time_plot(1,n);
            x_surf(:,(n-1)/10+1) = vertcat(x(:,n), xi_ox(1,n) + x_ox(:,n));
            xo_surf(:,(n-1)/10+1) = vertcat(xo(:,n), xo_ox(:,n));
        end
    end
    
    surf(time_surf, 1e4*x_surf, 100*xo_surf, 'EdgeColor', 'None');
    
%     semilogy(1e4*x(:, t_plot)', Co(:, t_plot)', 'r');
%     hold on;
%     plot([xi_ox(t_plot)*1e4 xi_ox(t_plot)*1e4], [Co_a_ox Co_ox_a], 'Color', 'r');
%     semilogy(1e4*(xi_ox(t_plot)+x_ox(:, t_plot))', Co_ox(:, t_plot)', 'r');
    
    colormap parula;
    xlim([0 time_plot(end)]);
    xlabel('Exposure time [s]');
    ylim([300 1e4*(e_Zr+d_total(end))]);
    ylabel('Cladding thickness [\mum]');
    zlim([0 100]);
    zlabel('Oxygen fration [%]')
    title('Evolution of oxygen fraction in the metal with exposure time');
end
%}


%% Oxygen fraction profile from the oxide-metal interface
%{

if strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')
    fig_xo = figure('Name', 'Oxygen fraction profile');
    
    % Load the EPMA oxygen profiles
    if strcmp(alloy, 'Zry4') && strcmp(temperature, '1100C')
        EPMA_file_list = {"100", "200", "500", "1000"};
        EPMA_color = {'r', [0 0.6 0], 'b', [0.5 0.5 0.5]};
    elseif strcmp(alloy, 'Zry4') && strcmp(temperature, '1200C')
        EPMA_file_list = {"50", "200", "500", "1000", "1500"};
        EPMA_color = {'r', '0 0.6 0', 'b', [0.5 0.5 0.5], [1 0.5 0]};
        
    end
    
    folder_path = '../../../Experiment_HT_corrosion_EPMA/EPMA/Blank_correction/';
    
    hold on;
%%{    
    for k = 1:length(EPMA_file_list)
        load(strcat(folder_path, 'EPMA_', temperature, '_', EPMA_file_list{k}, 's.mat'));
        t_plot = 10 * str2double(EPMA_file_list{k});
        
        O_both = vertcat(O_trav1, O_trav2);
        O_err_both = vertcat(O_err_trav1, O_err_trav2);
        [y_sort, idx_sort] = sort(vertcat(y_trav1, y_trav2));
        O_sort = O_both(idx_sort);
        O_err_sort = O_err_both(idx_sort);
        
        idx_xi_ox = find(y_sort > 0, 1);
        idx_xi_ab = find(y_sort > 1e4*(xi_ox(t_plot)-xi_ab(t_plot)), 1);
        
        O_smooth = vertcat(O_sort(1:idx_xi_ox-1,1), smooth(O_sort(idx_xi_ox:idx_xi_ab-1,1), 5, 'moving'), smooth(O_sort(idx_xi_ab:end,1), 5, 'moving'));

        errorbar(y_sort, O_smooth, O_err_sort, '-^', 'Color',  EPMA_color{k}, 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor', EPMA_color{k}, 'DisplayName', sprintf('t = %d s (EPMA)', 0.1*t_plot));
        
%         errorbar(y_trav1, O_trav1-5.5, O_err_trav1, '-^', 'Color',  EPMA_color{k}, 'MarkerFaceColor', EPMA_color{k}, 'DisplayName', sprintf('t = %d s (EPMA)', 0.1*t_plot));
%         errorbar(y_trav2, O_trav2-5.5, O_err_trav2, '-v', 'Color',  EPMA_color{k}, 'MarkerFaceColor',  EPMA_color{k}, 'DisplayName', sprintf('t = %d s (EPMA)', 0.1*t_plot));
        plot(1e4*(xi_ox(t_plot) - x(:, t_plot)), 100*xo(:, t_plot), 'Color',  EPMA_color{k}, 'DisplayName', sprintf('t = %d s (C4)', 0.1*t_plot));
%         plot([0 0], 100*[xo_ox(1, t_plot) xo_a_ox], 'Color',  EPMA_color{k}, 'DisplayName', sprintf('t = %d s (EPMA)', 0.1*t_plot));
%         plot(-1e4*x_ox(:, t_plot), 100*fliplr(xo_ox(:, t_plot)), 'Color',  EPMA_color{k}, 'DisplayName', sprintf('t = %d s (EPMA)', 0.1*t_plot));
        
    end
%%}
%{
    O_both = vertcat(O_trav1, O_trav2);
    O_err_both = vertcat(O_err_trav1, O_err_trav2);
    [y_sort, idx_sort] = sort(vertcat(y_trav1, y_trav2));
    O_sort = O_both(idx_sort);
    O_err_sort = O_err_both(idx_sort);
    
    load(strcat(folder_path, 'EPMA_', temperature, '_200s.mat'));
    errorbar(y_sort, O_sort, O_err_sort, '-^', 'Color', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 't = 200 s (EPMA)');
%}
    legend('-DynamicLegend', 'Location', 'NorthEast', 'AutoUpdate', 'Off');     % Legend is put before the last data to prevent them to appear in the legend (using AutoUpdate option)
    
    plot([0 0], 100*[xo_ox(1, t_plot) xo_a_ox], 'Color',  EPMA_color{k});
    plot(-1e4*x_ox(:, t_plot), 100*fliplr(xo_ox(:, t_plot)), 'Color',  EPMA_color{k});
    
    xlim([-20 200]);
    xlabel('Distance from the Oxide-Metal interface [\mum]');
    ylim([0 70]);
    ylabel('Oxygen fraction [%]');
%     title(sprintf('Oxygen concentration in the metal at T = %d C', T_avg - Tk));

%     print(fig_wg, 'Zry-4_xo_1100C', '-depsc');
end

%}


%% Species concentration in the oxide
%{
fig_C_ox = figure('Name', 'Species concentration in the oxide');

t_plot = 200;

semilogy(1e4*x_ox(:, t_plot)', Cv_ox(:, t_plot)', 'Color', [1 0.5 0], 'DisplayName', 'Vacancy concentration');
hold on;
semilogy(1e4*x_ox(:, t_plot)', Ce_ox(:, t_plot)', 'b', 'DisplayName', 'Electron concentration');
if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
    semilogy(1e4*x_ox(:, t_plot)', Ch_ox(:, t_plot)', 'Color', [0 0.6 0], 'DisplayName', 'Hydrogen concentration');
end
semilogy(1e4*x_ox(:, t_plot)', Co_ox(:, t_plot)', 'r', 'DisplayName', 'Oxygen concentration');

xlim([0 1e4*d_protect(t_plot)]);
xlabel('Oxide thickness [\mum]');
ylim([1e17 1e23]);
ylabel('Concentration [cm^{-3}]');
title(sprintf('Species concentration in the oxide at t = %d s', t_plot));
legend('-DynamicLegend', 'Location', 'SouthWest');
%}


%% Amount of Nb needed to compensate the space charges
%{
fig_Nb_needed = figure('Name', 'Amount of Nb needed to compensate the space charges');

t_plot = 1000;
load(strcat('output_file/', file_name_C4), 'Nb_needed');

hold on;
plot(100*x_ox(:, t_plot)/x_ox(end, t_plot), 100*Nb_needed(:, t_plot), 'b', 'DisplayName', 'C4 model');
if strcmp(model, 'C4') == false
    load(strcat('output_file/', file_name), 'Nb_needed');
    plot(100*x_ox(:, t_plot)/x_ox(end, t_plot), 100*Nb_needed(:, t_plot), 'r', 'DisplayName', strcat(model, ' model'));
end
plot([0 100], [0.5 0.5], '--k', 'DisplayName', 'Concentration of Nb in Zr-0.5Nb');

xlim([0 100]);
xlabel('Normalized oxide depth');
ylim([0 3]);
ylabel('Nb needed [wt %]');
title('Amount of oxidized Nb needed to compensate the space charges');
legend('-DynamicLegend', 'Location', 'NorthWest');
%}


%% Evolution of the interfaces positions with time

%{
fig_interfaces_positions = figure('Name', 'Interface position evolution');

area(time_plot, duct_b * 1e4, 'FaceColor', [0 1 0], 'FaceAlpha', 0.6, 'DisplayName', 'ductile \beta-Zr phase');
hold on;
fill([time_plot, fliplr(time_plot)], [xi_ab, fliplr(duct_b)] * 1e4, [0 0.6 0], 'FaceAlpha', 0.6, 'DisplayName', 'brittle \beta-Zr phase');
if max(T) > T_transf
    fill([time_plot, fliplr(time_plot)], [xi_ox, fliplr(xi_ab)] * 1e4, [0 0.3 0], 'FaceAlpha', 0.6, 'DisplayName', 'Oxygen-saturated \alpha-Zr phase');
end
fill([time_plot, fliplr(time_plot)], [xi_ox + d_protect, fliplr(xi_ox)] * 1e4, [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'DisplayName', 'Protective oxide ZrO_2');
if length(transition) > 1
    fill([time_plot, fliplr(time_plot)], [xi_ox + d_total, fliplr(xi_ox + d_protect)] * 1e4, 'black', 'FaceAlpha', 0.6, 'DisplayName', 'Non-protective oxide ZrO_2');
end

axis([0 time_plot(N+1) 300 700]);
% title('Interfaces positions evolution with exposure time');
lgnd = legend('-DynamicLegend', 'Location', 'SouthWest');

xlabel(strcat('Exposure time [', time_unit, ']'));
ylabel('Interfaces positions [\mum]');
%}



%% Gradients at the interfaces with right alpha/beta velocity

%%{
fig_grads = figure('Name', 'Gradients at the interfaces');

%subplot(2,1,1);
hold on;

plot(time_plot, grad_aox/Czr, 'Color', [1 0.5 0], 'DisplayName', 'Gradient at \alpha/oxide interface');
plot(time_plot, grad_ab/Czr, 'b', 'DisplayName', 'Gradient at \alpha/\beta interface');
plot(time_plot, grad_ba/Czr, 'Color', [0 0.6 0], 'DisplayName', 'Gradient at \beta/\alpha interface');

plot(time_MOOSE, grad_aox_MOOSE,'LineStyle', '--','Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, ICconst)');
plot(time_MOOSE, grad_ab_MOOSE,'LineStyle', '--', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, ICconst)');
plot(time_MOOSE, grad_ba_MOOSE,'LineStyle', '--','Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, ICconst)'); 


xlim([0 tmax]);
xticks(0:100:tmax);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylim([0 150]);
pbaspect([1.185 1 1]);
ylabel('Gradients of Co/Czr at the interface [/cm]');
title(strcat('Gradients at the interfaces evolution with exposure time, T=',temperature));
legend('-DynamicLegend', 'Location', 'NorthEast');
%}


%% Interfaces position Matlab vs MOOSE

%%{
fig_grads = figure('Name', 'Interfaces position evolution with time');

%subplot(2,1,2);
hold on;
plot(time_plot, xi_ox*1e4, 'Color', [1 0.5 0], 'DisplayName', '\alpha/oxide interface position');
plot(time_plot, xi_ab*1e4, 'b', 'DisplayName', '\alpha/\beta interface position');
plot(time_MOOSE, xi_ox_MOOSE, 'Color',[0.929 0.694 0.125],'LineStyle','--', 'DisplayName', 'MOOSE \alpha/oxide interface position');
plot(time_MOOSE, xi_ab_MOOSE, 'Color',[0.301 0.745 0.933],'LineStyle','--', 'DisplayName', 'MOOSE \alpha/\beta interface position');

xlim([0 tmax]);
xticks(0:100:tmax);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylim([350 600]);
yticks(350:50:600);
ylabel('Interfaces position [\mum]');
pbaspect([1.185 1 1]);
title(strcat('Interfaces position evolution with exposure time, T=',temperature));
legend('-DynamicLegend', 'Location', 'SouthWest');
%}

%% NRMSE on oxide and alpha layer thickness, MOOSE vs Matlab
%%{
fig_nrmse = figure('Name', 'NRMSE on oxide and alpha layer thickness evolution with time');
hold on;
plot(time_MOOSE, 100*d_diff_MOOSE, 'Color',[1 0.5 0], 'DisplayName', 'Oxide thickness NRMSE (Matlab vs MOOSE)');
plot(time_MOOSE, 100*alpha_diff_MOOSE, 'b', 'DisplayName', '\alpha layer thickness NRMSE (Matlab vs MOOSE)');
plot(time_MOOSE, 100*instant_d_diff, 'Color',[1 0.5 0],'LineStyle',':', 'DisplayName', 'Oxide thickness instantaneous error (Matlab vs MOOSE)');
plot(time_MOOSE, 100*instant_alpha_diff, 'b:', 'DisplayName', '\alpha layer thickness instantaneous error (Matlab vs MOOSE)');

xlim([0 tmax]);
xticks(0:100:tmax);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylim([0 20]);
yticks(0:100);
ylabel('NRMSE [%]');
pbaspect([1.185 1 1]);
title(strcat('Oxide and alpha layer thickness NRMSE evolution with time, T=',temperature));
legend('-DynamicLegend', 'Location', 'NorthWest');
%}

%% Alpha/beta interface velocity
%%{
fig_nrmse = figure('Name', 'Alpha/Beta interface velocity evolution with time');
hold on;
plot(time_plot, -1e4*v_ab, 'r', 'DisplayName', '|v_{\alpha/\beta}| (Matlab, from gradients)');
plot(time_plot(1:N), -1e4*v_ab_obs, 'k', 'DisplayName', '|v_{\alpha/\beta}| (Matlab,from interface position)');
plot(time_MOOSE, -1e4*v_ab_MOOSE, 'b', 'DisplayName', '|v_{\alpha/\beta}| (MOOSE)');

xlim([0 tmax]);
xticks(0:100:tmax);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylim([0 1.5]);
yticks(0:0.05:1.5);
ylabel('Velocity (\mum/s)');
pbaspect([1.185 1 1]);
title(strcat('\alpha/\beta interface velocity evolution with time, T=',temperature));
legend('-DynamicLegend', 'Location', 'NorthEast');
%}

%% Vacancy flux as a function of delta
%{
fig_Jv_delta = figure('Name', 'Vacancy flux as a function of oxide thickness');

hold on;
plot(10000*d_total, 10000*Jv, 'DisplayName', 'Vacancy flux calculated in Matlab');
plot(d_MOOSE, Jv_MOOSE,'Linestyle','None','Color','r','Marker','+', 'DisplayName', 'Vacancy flux calculated in MOOSE');

xlim([0 10000*d_total(end)]);
xlabel(strcat('Oxide thickness [\mum]'));
ylim([0 10000*max(Jv)]);
ylabel('Vacancy flux [/m^2/s]');
title('Vacancy flux as a function of oxide thickness');
legend('-DynamicLegend', 'Location', 'NorthEast');
%}
%% Save the figures in EPS format
% print(fig_name, 'name', '-depsc', '-painters', '-r300');