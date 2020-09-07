%% parabolic_rate_constant
%%Comparison of the parabolic rate constants for different models

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close;

%% Set default graphic properties

set(0, 'DefaultFigureColor', 'White');
set(0, 'DefaultAxesFontSize', 40);
set(0, 'DefaultAxesGridLineStyle', '-');
set(0, 'DefaultAxesXGrid', 'on');
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesXColor', 'Black');
set(0, 'DefaultAxesYColor', 'Black');
set(0, 'DefaultAxesZColor', 'Black');
set(0, 'DefaultAxesGridColor', 'Black');
set(0, 'DefaultAxeslineWidth', 1);
set(0, 'DefaultLineLineWidth', 4);
set(0, 'DefaultLineMarkerSize', 15);
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0.035 1 0.885]);
set(0, 'DefaultFigurePaperType', 'usletter');
set(0, 'DefaultFigurePaperOrientation', 'landscape');


%% Define the constants

R = 8.314;      % gas constant [J/K/mol]
RR = 1.987;      % gas constant [cal/K/mol]
Tk = 273.15;
Mo = 16;
MZr = 91.22;


%% Temperature ranges for each model
N = 100;
T1 = 1273;
T2 = 1773;
T = linspace(T1,T2,N+1);
T_plot = 1000./T;

XTick_label = {'1273' '1373' '1473' '1573' '1673' '1773'};


%% C4 model

T_C4_exp = [1273 1373 1473 1573 1673 1773];
T_C4_exp_plot = 1000./T_C4_exp;

T_C4 = linspace(1273,1773,N+1);
T_C4_plot = 1000./T_C4;

% wg_C4_exp = 1e-6 * [0.258152234183417 0.425566368193860 0.659661095818471 0.961975585166691 1.35218867044440 1.83614600749682].^2;      % Mallett-Perkins O diff coef
wg_C4_exp = 1e-6 * [0.240779494158052 0.401771132012664 0.627582994819798 0.923254707837994 1.37309093193377 1.84212328856162].^2;      % alpha-only optimization

wg_C4_fit = polyfit(T_C4_exp_plot, log(wg_C4_exp), 1);
wg_C4 = exp(wg_C4_fit(2)) * exp(1e3 * RR * wg_C4_fit(1) ./ (RR * T_C4));

fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
fprintf('Weight gain parabolic constant: %.2e * exp(%.0f / RT)\n', exp(wg_C4_fit(2)), 1e3*RR*wg_C4_fit(1));

% d_C4_exp = [0.000141584659027315 0.000222561017709353 0.000326039230996191 0.000467017231721446 0.000633602942777712 0.000830548518341550].^2;      % Mallett-Perkins O diff coef
d_C4_exp = [0.000151709330643642 0.000236101663724324 0.000345843043170161 0.000486700182941791 0.000621577265705473 0.000827344305329159].^2;      % alpha-only optimization

d_C4_fit = polyfit(T_C4_exp_plot, log(d_C4_exp), 1);
d_C4 = exp(d_C4_fit(2)) * exp(1e3 * RR * d_C4_fit(1) ./ (RR * T_C4));

fprintf('Oxide thickness parabolic constant: %.2e * exp(%.0f / RT)\n', exp(d_C4_fit(2)), 1e3*RR*d_C4_fit(1));

% alpha_C4_exp = [0.000240809119231403 0.000329186069122622 0.000425806228881063 0.000703602103022664 0.000980770771809756 0.00133115451678848].^2;       % Mallett-Perkins O diff coef
alpha_C4_exp = [9.14827932026206e-05 0.000193377985907255 0.000371754268339528 0.000557328782614102 0.000889159297530432 0.00130492760838661].^2;       % alpha-only optimization

alpha_C4_fit = polyfit(T_C4_exp_plot, log(alpha_C4_exp), 1);
alpha_C4 = exp(alpha_C4_fit(2)) * exp(1e3 * RR * alpha_C4_fit(1) ./ (RR * T_C4));

fprintf('Alpha thickness parabolic constant: %.2e * exp(%.0f / RT)\n', exp(alpha_C4_fit(2)), 1e3*RR*alpha_C4_fit(1));


%% Cathcart-Pawel model

T_CP = linspace(1273,1773,N+1);
T_CP_plot = 1000./T_CP;

wg_CP = 2 * 0.1811 * exp(-39940./(RR*T_CP));

d_CP = 2 * 0.01126 * exp(-35890./(RR*T_CP));

alpha_CP = 2 * 0.7615 * exp(-48140./(RR*T_CP));


%% Baker-Just model

T_BJ = linspace(1273,2123,N+1);
T_BJ_plot = 1000./T_BJ;

wg_BJ = (2 * Mo / MZr)^2 * 33.3 * exp(-45500./(RR*T_BJ));

d_BJ = (2 * Mo / MZr)^2 / 1.477^2 * 33.3 * exp(-45500./(RR*T_BJ));


%% Lemmon experiment

T_Lemmon = linspace(1273,1963,N+1);
T_Lemmon_plot = 1000./T_Lemmon;

wg_Lemmon = (Mo / 22.71)^2 * 0.1132 * exp(-34000./(RR*T_Lemmon));


%% Biederman model

T_WPI = linspace(1144,1755,N+1);
T_WPI_plot = 1000./T_WPI;

wg_WPI = 0.1953^2 * exp(-33370./(RR*T_WPI));

d_WPI = 364.4e-4^2 * exp(-27840./(RR*T_WPI));

alpha_WPI = 2257.5e-4^2 * exp(-37690./(RR*T_WPI));


%% Kawasaki model

T_Kawasaki = linspace(1273,1603,N+1);
T_Kawasaki_plot = 1000./T_Kawasaki;

wg_Kawasaki = 0.468 * exp(-40710./(RR*T_Kawasaki));

d_Kawasaki = 0.0215 * exp(-35860./(RR*T_Kawasaki));

alpha_Kawasaki = (sqrt(0.306 * exp(-38970./(RR*T_Kawasaki))) - sqrt(d_Kawasaki)).^2;
alpha_Kawasaki_fit = 0.2183 * exp(-41508 ./ (RR*T_Kawasaki));


%% Leistikow model

T_Leistikow = linspace(1073,1773,N+1);
T_Leistikow_plot = 1000./T_Leistikow;

wg_Leistikow = 5.24e-1 * exp(-174300 ./ (R*T_Leistikow));

d_Leistikow = 7.82e-2 * exp(-168100./(R*T_Leistikow));


%% Urbanic model

T_Urbanic = linspace(1323,1853,N+1);
T_Urbanic_plot = 1000./T_Urbanic;

wg_Urbanic = (2 * Mo / MZr)^2 * 0.296 * exp(-16820./T_Urbanic);

d_Urbanic = 0.036^2 * exp(-27000./(RR*T_Urbanic));

alpha_Urbanic = 0.390^2 * exp(-39400./(RR*T_Urbanic));


%% Others

% wg_Hobson = 0.1553 * exp(-39290./(RR*T));
% wg_Klepfer = 0.02203 * exp(-33500./(RR*T));


%% Plot weight gain

fig_wg = figure();

semilogy(T_C4_exp_plot, wg_C4_exp, '^r', 'MarkerFaceColor', 'r', 'DisplayName', 'C4 model');
hold on;
semilogy(T_CP_plot, wg_CP, 'k', 'DisplayName', 'Cathcart-Pawel model: 0.362 * exp(-39940 / RT)');
semilogy(T_BJ_plot, wg_BJ, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model: 4.10 * exp(-45500 / RT)');
semilogy(T_WPI_plot, wg_WPI, 'g', 'DisplayName', 'Biederman model: 3.81e-02 * exp(-33370 / RT)');
semilogy(T_Kawasaki_plot, wg_Kawasaki, 'Color', [0 0.6 0], 'DisplayName', 'Kawasaki model: 0.468 * exp(-40710 / RT)');
semilogy(T_Leistikow_plot, wg_Leistikow, 'Color', [0.5 0 0.5], 'DisplayName', 'Leistikow model: 0.524 * exp(-41660 / RT)');
semilogy(T_Urbanic_plot, wg_Urbanic, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Urbanic model: 3.64e-02 * exp(-33421 / RT)');
% semilogy(T_Lemmon_plot, wg_Lemmon, 'm', 'DisplayName', 'Lemmon model');
% semilogy(T_plot, wg_Hobson, 'c', 'DisplayName', 'Hobson model');
% semilogy(T_plot, wg_Klepfer, 'b', 'DisplayName', 'Klepfer model');

semilogy(T_C4_plot, wg_C4, 'r', 'DisplayName', sprintf('C4 model: %.2e * exp(%.0f / RT)\n', exp(wg_C4_fit(2)), 1e3*RR*wg_C4_fit(1)));

% Create the bottom axis
ax1 = gca;
xlim([1000/T2 1000/T1]);
% ylim([1e-9 1e-4]);
xlabel('1000/T [K^{-1}]');
ylabel('Weight gain parabolic rate constant [(g/cm^2)^2/s]');
legend('-DynamicLegend', 'Location', 'SouthWest');
grid off;

% Create the top axis
XTick_vect = linspace(T1,T2,length(XTick_label));
for k = 2:length(XTick_vect)-1
    XTick_vect(k) = T1+(1/T1 - 1/XTick_vect(k))/(1/T1 - 1/T2)*(T2-T1);
end
set(ax1, 'YGrid', 'On', 'YMinorGrid', 'On', 'Box', 'off');         % suppress the box surrounding the plot
ax1_pos = get(ax1, 'Position');     % extract the position
ax1_pos(2) = 0.82 * ax1_pos(2);
ax1_pos(4) = 1.03 * ax1_pos(4);
set(ax1, 'Position', ax1_pos);
ax2 = axes('Position', ax1_pos);    % create the new axis
xlim([T(1) T(N+1)]);
set(ax2, 'XDir', 'Reverse', 'Color', 'None', 'YTick', [], 'XGrid', 'On', 'YGrid', 'Off', 'Box', 'Off', 'XAxisLocation', 'Top', 'XTick', XTick_vect, 'XTickLabel', XTick_label);
xlabel('T [K]')

% Add a fake axis on the right to close the box
axes('Color', 'none', 'XTick', [], 'YTick', [], 'XGrid', 'Off', 'YGrid', 'Off', 'YAxisLocation', 'Right');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot oxide thickness

fig_d = figure();

semilogy(T_C4_exp_plot, d_C4_exp, '^r', 'MarkerFaceColor', 'r', 'DisplayName', 'C4 model');
hold on;
semilogy(T_CP_plot, d_CP, '-k', 'DisplayName', 'Cathcart-Pawel model: 2.25e-02 * exp(-35890 / RT)');
semilogy(T_BJ_plot, d_BJ, 'Color', [1 0.5 0], 'DisplayName', 'Baker-Just model: 1.88 * exp(-45500 / RT)');
semilogy(T_WPI_plot, d_WPI, '-g', 'DisplayName', 'Biederman model: 1.33e-03 * exp(-27840 / RT))');
semilogy(T_Kawasaki_plot, d_Kawasaki, 'Color', [0 0.6 0], 'DisplayName', 'Kawasaki model: 2.15e-02 * exp(-35860 / RT)');
semilogy(T_Leistikow_plot, d_Leistikow, 'Color', [0.5 0 0.5], 'DisplayName', 'Leistikow model: 7.82e-02 * exp(-40175 / RT)');
semilogy(T_Urbanic_plot, d_Urbanic, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Urbanic model: 1.30e-03 * exp(-27000 / RT)');
semilogy(T_C4_plot, d_C4, '-r', 'DisplayName', sprintf('C4 model: %.2e * exp(%.0f / RT)\n', exp(d_C4_fit(2)), 1e3*RR*d_C4_fit(1)));

% Create the bottom axis
ax1 = gca;
xlim([1000/T2 1000/T1]);
% ylim([1e-8 1e-5]);
xlabel('1000/T [K^{-1}]');
ylabel('Oxide thickness parabolic rate constant [(cm)^2/s]');
legend('-DynamicLegend', 'Location', 'NorthEast');
grid off;

% Create the top axis
XTick_vect = linspace(T1,T2,length(XTick_label));
for k = 2:length(XTick_vect)-1
    XTick_vect(k) = T1+(1/T1 - 1/XTick_vect(k))/(1/T1 - 1/T2)*(T2-T1);
end
set(ax1, 'YGrid', 'On', 'YMinorGrid', 'On', 'Box', 'off');         % suppress the box surrounding the plot
ax1_pos = get(ax1, 'Position');     % extract the position
ax1_pos(2) = 0.82 * ax1_pos(2);
ax1_pos(4) = 1.03 * ax1_pos(4);
set(ax1, 'Position', ax1_pos);
ax2 = axes('Position', ax1_pos);    % create the new axis
xlim([T(1) T(N+1)]);
set(ax2, 'XDir', 'Reverse', 'Color', 'None', 'YTick', [], 'XGrid', 'On', 'YGrid', 'Off', 'Box', 'Off', 'XAxisLocation', 'Top', 'XTick', XTick_vect, 'XTickLabel', XTick_label);
xlabel('T [K]')

% Add a fake axis on the right to close the box
axes('Color', 'none', 'XTick', [], 'YTick', [], 'XGrid', 'Off', 'YGrid', 'Off', 'YAxisLocation', 'Right');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot alpha thickness

fig_alpha = figure();

semilogy(T_C4_exp_plot, alpha_C4_exp, '^r', 'MarkerFaceColor', 'r', 'DisplayName', 'C4 model');
hold on;
semilogy(T_CP_plot, alpha_CP, '-k', 'DisplayName', 'Cathcart-Pawel model: 1.52 * exp(-48140 / RT)');
semilogy(T_WPI_plot, alpha_WPI, '-g', 'DisplayName', 'Biederman model: 5.10e-02 * exp(-37690 / RT)');
semilogy(T_Kawasaki_plot, alpha_Kawasaki, 'Color', [0 0.6 0], 'DisplayName', 'Kawasaki model: 0.218 * exp(-41508 / RT)');
semilogy(T_Urbanic_plot, alpha_Urbanic, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Urbanic model: 0.152 * exp(-39400 / RT)');
semilogy(T_C4_plot, alpha_C4, '-r', 'DisplayName', sprintf('C4 model: %.2f * exp(%.0f / RT)\n', exp(alpha_C4_fit(2)), 1e3*RR*alpha_C4_fit(1)));

% Create the bottom axis
ax1 = gca;
xlim([1000/T2 1000/T1]);
% ylim([1e-8 1e-5]);
xlabel('1000/T [K^{-1}]');
ylabel('\alpha thickness parabolic rate constant [(cm)^2/s]');
legend('-DynamicLegend', 'Location', 'SouthWest');
grid off;

% Create the top axis
XTick_vect = linspace(T1,T2,length(XTick_label));
for k = 2:length(XTick_vect)-1
    XTick_vect(k) = T1+(1/T1 - 1/XTick_vect(k))/(1/T1 - 1/T2)*(T2-T1);
end
set(ax1, 'YGrid', 'On', 'YMinorGrid', 'On', 'Box', 'off');         % suppress the box surrounding the plot
ax1_pos = get(ax1, 'Position');     % extract the position
ax1_pos(2) = 0.82 * ax1_pos(2);
ax1_pos(4) = 1.03 * ax1_pos(4);
set(ax1, 'Position', ax1_pos);
ax2 = axes('Position', ax1_pos);    % create the new axis
xlim([T(1) T(N+1)]);
set(ax2, 'Position', ax1_pos, 'XDir', 'Reverse', 'Color', 'None', 'YTick', [], 'XGrid', 'On', 'YGrid', 'Off', 'Box', 'Off', 'XAxisLocation', 'Top', 'XTick', XTick_vect, 'XTickLabel', XTick_label);
xlabel('T [K]')

% Add a fake axis on the right to close the box
axes('Color', 'none', 'XTick', [], 'YTick', [], 'XGrid', 'Off', 'YGrid', 'Off', 'YAxisLocation', 'Right');

%% Save the figures in EPS format
% print(fig_wg, 'empirical_model_wg', '-depsc');
% print(fig_d, 'empirical_model_d', '-depsc');
% print(fig_alpha, 'empirical_model_alpha', '-depsc');