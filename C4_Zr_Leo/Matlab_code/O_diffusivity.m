%% O_diffusivity
%%Comparison of the oxygen diffusion coefficient in literature

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set default graphic properties

set(0, 'DefaultFigureColor', 'White');
set(0, 'DefaultAxesFontSize', 30);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define the constants

R = 8.314;      % gas constant [J/K/mol]
RR = 1.987;      % gas constant [cal/K/mol]
Tk = 273;


%% Temperature range
N = 100;
T1 = 1273;
T2 = 1773;
T = linspace(T1,T2,N+1);
T_plot = 1000./T;

XTick_label = {'1273' '1373' '1473' '1573' '1673' '1773'};
% XTick_label = {'673' '773' '873' '973' '1073' '1173' '1273' '1373' '1473' '1573' '1673' '1773'};      % Use if temperature range 673-1773K

%% Zanella

D_Zanella = exp(0.6282) * exp(-24910./T);


%% Ritchie (290C-1500C)

T_Ritchie = linspace(923,1773,N+1);
T_Ritchie_plot = 1000./T_Ritchie;

% Low temperature : 290C < T < 650C
Da_Ritchie_LT = 0.0661 * exp(-44000./(RR*T_Ritchie));

% High temperature : 650C < T < 1500C
Da_Ritchie_HT = 16.5 * exp(-54700./(RR*T_Ritchie));

Da_Ritchie = zeros(1,N+1);
for n = 1:N+1
    if T(n) < 650 + 273
        Da_Ritchie(n) = Da_Ritchie_LT(n);
    else
        Da_Ritchie(n) = Da_Ritchie_HT(n);
    end
end


%% Perkins (900C-1500C) == Pawel

T_Perkins = linspace(1173,1773,N+1);
T_Perkins_plot = 1000./T_Perkins;
T_Perkins_plot_uplow = [T_Perkins_plot, fliplr(T_Perkins_plot)];

Db_Perkins = 2.48e-2 * exp(-28200./(RR*T_Perkins));

Db_Perkins_up = 3.60e-2 * exp(-27100./(RR*T_Perkins));
Db_Perkins_low = 1.71e-2 * exp(-29300./(RR*T_Perkins));

Db_Perkins_uplow = [Db_Perkins_up, fliplr(Db_Perkins_low)];

% Experimental
T_Perkins_exp = [905 930 955 952 1001 1005 1048 1049 1102 1105 1150 1153 1196 1198 1199 1200 1248 1250 1300 1309 1346 1347 1397 1398 1445 1450 1500 1504 1502 ...
                 1099 1102 1149 1152 1199 1200 1201 1205 1249 1250 1296 1301 1347 1348 1398 1401 1446 1449] + 273;
T_Perkins_exp_plot = 1000./T_Perkins_exp;
Db_Perkins_exp = [0.254 0.674 1.96 1.50 3.01 3.76 5.90 5.95 6.40 7.53 12.8 10.1 19.1 16.4 13.6 14.4 26.0 25.2 21.2 30.1 35.3 50.0 47.7 58.8 66.9 73.0 92.0 71.4 81.3 ...
                  9.35 8.01 12.1 10.4 19.3 20.2 16.5 18.3 22.8 22.3 30.0 29.4 40.6 36.1 48.1 51.0 60.6 58.0] * 1e-7;


%% Mallett (1000C-1500C)

T_Mallett = linspace(1273,1773,N+1);
T_Mallett_plot = 1000./T_Mallett;
T_Mallett_plot_uplow = [T_Mallett_plot, fliplr(T_Mallett_plot)];

Da_Mallett = 0.196 * exp(-41000./(RR*T_Mallett));
Da_Mallett_low = 0.196 * exp(-42500./(RR*T_Mallett));
Da_Mallett_up = 0.196 * exp(-39500./(RR*T_Mallett));
Da_Mallett_uplow = [Da_Mallett_low, fliplr(Da_Mallett_up)];

Db_Mallett = 0.0453 * exp(-28200./(RR*T_Mallett));
Db_Mallett_low = 0.0453 * exp(-30600./(RR*T_Mallett));
Db_Mallett_up = 0.0453 * exp(-25800./(RR*T_Mallett));
Db_Mallett_uplow = [Db_Mallett_low, fliplr(Db_Mallett_up)];

% Experimental
T_Mallett_exp = [1000 1100 1200 1300 1400 1500] + 273;
T_Mallett_exp_plot = 1000./T_Mallett_exp;
Da_Mallett_exp = [1.8e-8 4.8e-8 2.3e-7 3.6e-7 7.9e-7 1.7e-6];
Db_Mallett_exp = [8.2e-7 1.4e-6 2.4e-6 3.8e-6 8.1e-6 2.4e-5];


%% Pawel (1000C-1500C) == Perkins

T_Pawel = linspace(1273,1773,N+1);
T_Pawel_plot = 1000./T_Pawel;
T_Pawel_plot_uplow = [T_Pawel_plot, fliplr(T_Pawel_plot)];

Db_Pawel = 0.0263 * exp(-28200./(RR*T_Pawel));
Db_Pawel_up = 1.45*0.0263 * exp(-0.96*28200./(RR*T_Pawel));
Db_Pawel_low = 0.69*0.0263 * exp(-1.04*28200./(RR*T_Pawel));
Db_Pawel_uplow = [Db_Pawel_low, fliplr(Db_Pawel_up)];


%% C4 model

T_C4_exp = 1000 ./ [1273 1373 1473 1573 1673 1773];
% Alpha only optimization
Da_C4_exp = [4.5e-09 2.6e-08 1e-7 2.61152343691429e-07 8.72249840061976e-07 1.73131961527490e-06];
Db_C4_exp = [3.8e-07 3.7e-07 6e-7 2.01070471839367e-06 8.56187744202073e-06 9.27853860803972e-06];

T_C4 = linspace(673,1773,N+1);
T_C4_plot = 1000./T_C4;

Da_C4_fit = polyfit(T_C4_exp, log(Da_C4_exp), 1);
Db_C4_fit = polyfit(T_C4_exp, log(Db_C4_exp), 1);
Da_C4 = exp(Da_C4_fit(2)) * exp(1e3 * RR * Da_C4_fit(1) ./ (RR * T_C4));
Db_C4 = exp(Db_C4_fit(2)) * exp(1e3 * RR * Db_C4_fit(1) ./ (RR * T_C4));


%% Corvolan

T_Corvolan_exp = [1373 1473 1523];
T_Corvolan_exp_plot = 1000./T_Corvolan_exp;
Da_Corvolan_exp = [2.48e-8 8.50e-8 1.48e-7];
Db_Corvolan_exp = [8.44e-7 1.55e-6 2.05e-6];

T_Corvolan = linspace(1373,1523,N+1);
T_Corvolan_plot = 1000./T_Corvolan;

Da_Corvolan = 1.88 * exp(-49500 ./ (RR*T_Corvolan));
Db_Corvolan = 6.537e-3 * exp(-24430 ./ (RR*T_Corvolan));


%% Pemsler (400C-1500C)

T_Pemsler = linspace(673,1773,N+1);
T_Pemsler_plot = 1000./T_Pemsler;

% Low temperature: 400C < T < 585C
Da_Pemsler_LT = 9.4 * exp(-51780./(RR*T_Pemsler));

% Least square with high temperature data from Mallett: 400C < T < 1500C
Da_Pemsler = 5.2 * exp(-50800./(RR*T_Pemsler));


%% Zhang (300C-450C & 175C-275C)

% Low temperature on alpha Zr-1.0Nb: 300C < 450C
Da_Zhang = 0.172 * exp(-187.47e3./(R*T));

% Low temperature on beta Zr-20Nb: 175C < T < 275C
Db_Zhang = 0.69 * exp(-149.45e3./(R*T));


%% plot the figure

fig = figure;

% Da
semilogy(T_Mallett_exp_plot, Da_Mallett_exp * 1e-4, '^g', 'MarkerFaceColor', 'g');
hold on;
alpha_Mallett = semilogy(T_Mallett_plot, Da_Mallett * 1e-4, ':g', 'DisplayName', 'D_\alpha Mallett: 0.196 * exp(-41000 / RT)');
fill(T_Mallett_plot_uplow, Da_Mallett_uplow * 1e-4, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'w');

% alpha_Ritchie = semilogy(T_Ritchie_plot, Da_Ritchie * 1e-4, 'Color', [0 0.6 0], 'DisplayName', 'Ritchie');
alpha_Pemsler = semilogy(T_Pemsler_plot, Da_Pemsler * 1e-4, 'Color', [0.5 0.5 0.5], 'DisplayName', 'D_\alpha Pemsler: 5.2 * exp(-50800 / RT)');

semilogy(T_Corvolan_exp_plot, Da_Corvolan_exp * 1e-4, '^b', 'MarkerFaceColor', 'b');
alpha_Corvolan = semilogy(T_Corvolan_plot, Da_Corvolan * 1e-4, ':b', 'DisplayName', 'D_\alpha Corvolan: 1.88 * exp(-49500 / RT)');

semilogy(T_C4_exp, Da_C4_exp * 1e-4, '^r', 'MarkerFaceColor', 'r');
alpha_C4 = semilogy(T_C4_plot, Da_C4 * 1e-4, ':r', 'DisplayName', sprintf('D_{\\alpha} C4 model: %.2f * exp(%.0f / RT)\n', exp(Da_C4_fit(2)), 1e3*RR*Da_C4_fit(1)));

% Db
semilogy(T_Mallett_exp_plot, Db_Mallett_exp * 1e-4, 'vg', 'MarkerFaceColor', 'g');
beta_Mallett = semilogy(T_Mallett_plot, Db_Mallett * 1e-4, 'g', 'DisplayName', 'D_\beta Mallett: 4.53e-02 * exp(-28200 / RT)');
fill(T_Mallett_plot_uplow, Db_Mallett_uplow * 1e-4, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'w');

semilogy(T_Perkins_exp_plot, Db_Perkins_exp * 1e-4, 'vm', 'MarkerFaceColor', 'm');
beta_Perkins = semilogy(T_Perkins_plot, Db_Perkins * 1e-4, 'm', 'DisplayName', 'D_\beta Perkins: 2.48e-02 * exp(-28200 / RT)');
fill(T_Perkins_plot_uplow, Db_Perkins_uplow * 1e-4, 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'w');

semilogy(T_Corvolan_exp_plot, Db_Corvolan_exp * 1e-4, 'vb', 'MarkerFaceColor', 'b');
beta_Corvolan = semilogy(T_Corvolan_plot, Db_Corvolan * 1e-4, 'b', 'DisplayName', 'D_\beta Corvolan: 6.54e-03 * exp(-24430 / RT)');

semilogy(T_C4_exp, Db_C4_exp * 1e-4, 'vr', 'MarkerFaceColor', 'r');
beta_C4 = semilogy(T_C4_plot, Db_C4 * 1e-4, 'r', 'DisplayName', sprintf('D_{\\beta} C4 model: %.2e * exp(%.0f / RT)\n', exp(Db_C4_fit(2)), 1e3*RR*Db_C4_fit(1)));


% Create the bottom axis
ax1 = gca;
xlim([1000/T2 1000/T1]);
xlabel('1000/T [K^{-1}]');
ylim([1e-13 1e-8]);
ylabel('Oxygen diffusion coefficient [m^2/s]');
legend([alpha_C4, alpha_Corvolan, alpha_Mallett, alpha_Pemsler, beta_C4, beta_Corvolan, beta_Mallett, beta_Perkins], 'Location', 'SouthWest');
grid off;

% Create the top axis
XTick_vect = linspace(T1,T2,length(XTick_label));
for k = 2:length(XTick_vect)-1
    XTick_vect(k) = T1+(1/T1 - 1/XTick_vect(k))/(1/T1 - 1/T2)*(T2-T1);
end
set(ax1, 'YGrid', 'On', 'YMinorGrid', 'On', 'Box', 'off');         % suppress the box surrounding the plot
ax1_pos = get(ax1, 'Position');     % extract the position
% ax1_pos(2) = 0.82 * ax1_pos(2);
% ax1_pos(4) = 1.03 * ax1_pos(4);
set(ax1, 'Position', ax1_pos);
ax2 = axes('Position', ax1_pos);    % create the new axis
xlim([T(1) T(N+1)]);
set(ax2, 'Position', ax1_pos, 'XDir', 'Reverse', 'Color', 'None', 'YTick', [], 'XGrid', 'On', 'YGrid', 'Off', 'Box', 'Off', 'XAxisLocation', 'Top', 'XTick', XTick_vect, 'XTickLabel', XTick_label);
xlabel('T [K]')

% Add a fake axis on the right to close the box
axes('Color', 'none', 'XTick', [], 'YTick', [], 'XGrid', 'Off', 'YGrid', 'Off', 'YAxisLocation', 'Right');

%% Save the figures in EPS format
print(fig, 'O_diff_coef', '-depsc');