%% C4_exp_data
%%Implementation of the different experimental values and their fit

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set default graphic properties

set(0, 'DefaultFigureColor', 'White');
set(0, 'DefaultAxesFontSize', 30);      % 50 for presentation ; 30 else
set(0, 'DefaultAxesGridLineStyle', '-');
set(0, 'DefaultAxesXGrid', 'on');
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesXColor', 'Black');
set(0, 'DefaultAxesYColor', 'Black');
set(0, 'DefaultAxesZColor', 'Black');
set(0, 'DefaultAxeslineWidth', 1);
set(0, 'DefaultLineLineWidth', 3);      % 5 for presentation ; 3 else
set(0, 'DefaultLineMarkerSize', 10);
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0.03 1 0.87]);
set(0, 'DefaultFigurePaperType', 'usletter');
set(0, 'DefaultFigurePaperOrientation', 'landscape');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Operating Temperature experimental data

%% Experimental values for Zr-0.5Nb with C4_360C_240d

if strcmp(alloy, 'Zr05Nb') && strcmp(temperature, '360C')
	day_wg_exp = [3, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240];%, 255, 265];
	wg_exp =  [7.21, 14.38, 19.22, 23.62, 26.62, 30.25, 32.11, 35.17, 36.75, 39.39, 40.72, 43.55, 44.69, 47.42, 48.78, 52.68, 52.57];%, 56.88, 60.29];
	d_exp = wg_exp./14.77;

	wg_fit = 3.735 * days.^0.4815;

	day_Ch_exp = [0, 15, 30, 45, 60, 75, 90, 105, 135, 165, 195, 240];
	Ch_exp = [14.37, 15.51, 15.22, 16.81, 16.29, 17.15, 17.54, 16.66, 21.70, 27.28, 38.85, 61.25];
    
    days_fit = (days - 96.25) ./ 74.72;
	Ch_fit = -0.5439 * days_fit.^5 + 0.1066 * days_fit.^4 + 4.525 * days_fit.^3 + 5.066 * days_fit.^2 + 3.25 * days_fit + 17.35;

	day_fh_exp = [0, 15, 30, 45, 60, 75, 90, 105, 135, 165, 195, 240];
	fh_exp = [0, 0.00935, 0.00494, 0.01158, 0.00852, 0.01077, 0.01136, 0.00770, 0.02162, 0.03362, 0.06073, 0.10409];
    fh_exp_err = [0, 0.00008, 0.00017, 0.00041, 0.00057, 0.00108, 0.00033, 0.00015, 0.00078, 0.00085, 0.00071, 0];

	fh_fit = -9.299e-11 * days.^4 + 5.464e-8 * days.^3 - 7.567e-6 * days.^2 + 0.0003747 * days + 0.003375;

%% Experimental values for Zr-1.0Nb with C4_360C_225d or C4_360C_60d

elseif strcmp(alloy, 'Zr10Nb') && strcmp(temperature, '360C')
% 	day_wg_exp = [1.00, 3.00, 7.00, 15.00, 15.00, 15.00, 15.00, 30.00, 30.00, 30.00, 30.00, 45.00, 45.00, 45.00, 45.00, 60.00, 60.00, 60.00, 72.00, 75.00, 75.00, 75.00, 87.00, 90.00, 90.00, 102.00, 105.00, 105.00, 117.00, 120.00, 120.00, 135.00 ,135.00, 150.00, 150.00, 165.00, 165.00, 180.00, 180.00, 195.00, 195.00, 210.00, 210.00, 225.00];
% 	wg_exp =  [6.25, 10.01, 13.95, 20.23, 20.74, 21.11, 21.02, 28.36, 28.73, 29.37, 29.10, 35.12, 35.31, 36.16, 35.75, 42.48, 42.23, 42.48, 49.22, 52.53, 51.99, 52.17, 63.79, 66.64, 66.57, 74.11, 75.74, 75.82, 82.49, 84.04, 83.19, 91.67, 90.10, 105.65, 102.29, 114.39, 111.67, 123.38, 120.66, 128.08, 128.20, 136.70, 138.28, 150.61]; 
	
    day_wg_exp = cat(2, 1*ones(1,8), 3*ones(1,7), 7*ones(1,6), 15*ones(1,15), 30*ones(1,14), 45*ones(1,13), ...
        60*ones(1,10), 72*ones(1,2), 75*ones(1,8), 87*ones(1,3), 90*ones(1,4), 102*ones(1,4), 105*ones(1,6), ...
        117, 120*ones(1,3), 135*ones(1,5), 150*ones(1,2), 165*ones(1,4), 180*ones(1,3), 195*ones(1,2), 210*ones(1,2), 225);
    wg_exp = [06.4, 06.4, 06.1, 06.3, 06.1, 06.1, 06.1, 06.5, ... 1
              10.2, 10.1, 10.1, 09.7, 09.7, 10.0, 10.2 ... 3
              13.7, 14.3, 13.8, 13.7, 14.0, 14.3, ... 7
              20.8, 19.7, 19.7, 20.1, 20.9, 19.9, 21.4, 20.8, 21.4, 20.1, 20.7, 21.2, 21.5, 21.1, 21.0, ... 15
              27.9, 27.8, 28.3, 29.4, 27.8, 29.6, 29.2, 29.5, 27.6, 28.9, 29.4, 29.8, 29.3, 28.9, ... 30
              34.3, 35.0, 36.1, 34.3, 36.4, 35.7, 36.2, 33.9, 35.8, 36.0, 36.7, 35.8, 35.7, ... 45
              41.8, 43.2, 41.5, 42.9, 42., 42.8, 41.5, 41.9, 42.4, 43.1, ... 60
              49.3, 49.1, ... 72
              52.5, 52.3, 51.4, 52.2, 51.3, 52.6, 52.0, 52.4, ... 75
              63.6, 64.3, 63.5, ... 87
              66.6, 66.4, 66.7, 66.6, ... 90
              74.6, 74.4, 74.5, 72.9, ... 102
              76.1, 76.0, 76.3, 74.5, 75.6, 76.0, ... 105
              82.5, ... 117
              84.0, 82.9, 83.4, ... 120
              92.2, 90.5, 92.3, 90.3, 89.9, ... 135
              105.7, 102.3, ... 150
              115.0, 113.7, 111.0, 112.3, ... 165
              123.4, 120.8, 120.5, ... 180
              128.1, 128.2, ... 195
              136.7, 138.3, ... 210
              150.6, ... 225
              ];
          
    d_exp = wg_exp/14.77;
    
	wg_fit = 5.256 * days.^0.5055;

	day_Ch_exp = [1, 1, 3, 3, 7, 7, 15, 15, 30, 45, 45, 60, 60, 75, 75, 90, 90, 105, 105, 120, 120, 135, 135, 150, 150, 165, 165, 180, 195, 210, 225];

	Ch_exp = [7.2, 7.25, 9.26, 9.76, 10.96, 10.31, 12.27, 11.56, 17.93, 27.16, 28.44, 40.43, 40.94, 54.26, 53.34, 69.77, 69.70, 79.6, 79.75, 90.68, 87.46, 101.11, 101.41, 116.84, 115.66, 126.26, 123.85, 132.79, 150.60, 162.69, 184.94];

    days_fit = (days - 89.74) ./ 68.9;
	Ch_fit = -2.134 * days_fit.^8 + 7.549 * days_fit.^7 + 0.5857 * days_fit.^6 - 22.76 * days_fit.^5 + 14.03 * days_fit.^4 + 13.63 * days_fit.^3 - 12.89 * days_fit.^2 + 55.01 * days_fit + 67.57;

	day_fh_exp = [1, 3, 7, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225];
	fh_exp = [0, 0.01333, 0.01579, 0.01493, 0.02580, 0.04167, 0.05655, 0.06225, 0.06711, 0.06867, 0.07072, 0.07325, 0.07744, 0.07451, 0.07489, 0.08174, 0.08122, 0.08525];
%     fh_exp = [0.079, 0.065, 0.055, 0.041, 0.045, 0.057, 0.069, 0.072, 0.075, 0.076, 0.077, 0.079, 0.083, 0.079, 0.079, 0.086, 0.085, 0.089];
    fh_exp_err = [0.0057, 0.0042, 0.0034, 0.0024, 0.0024, 0.0030, 0.0035, 0.0036, 0.0038, 0.0038, 0.0038, 0.0039, 0.0041, 0.0039, 0.0039, 0.0042, 0.0042, 0.0044];

	fh_fit = 7.948e-15 * days.^6 - 6.962e-12 * days.^5 + 2.335e-9 * days.^4 - 3.574e-7 * days.^3 + 2.098e-5*days.^2 + 0.000383*days + 0.006872;

end


%% IAEA-CRP
 
if strcmp(alloy, 'Zr10Nb')
    if strcmp(temperature, '500C')
        time_exp_IAEA_CRP = [163000 653000 653000] / 86400;
        wg_exp_IAEA_CRP = [0.461 0.986 0.980];
        k_exp_IAEA_CRP = [0.00114 0.00122 0.00121];
    elseif strcmp(temperature, '600C')
        time_exp_IAEA_CRP = [2400 2400 10800 138600] / 86400;
        wg_exp_IAEA_CRP = [0.276 0.435 0.502 1.980];
        k_exp_IAEA_CRP = [0.00563 0.00888 0.00483 0.00532];
    elseif strcmp(temperature, '700C')
        time_exp_IAEA_CRP = [171 171 1200 4800 4800];
        wg_exp_IAEA_CRP = [0.245 0.261 0.635 1.060 1.058];
        k_exp_IAEA_CRP = [0.01875 0.01994 0.01834 0.01530 0.01528];
    elseif strcmp(temperature, '800C')
        time_exp_IAEA_CRP = [69 69 277 277 1107 1107];
        wg_exp_IAEA_CRP = [0.390 0.380 0.796 0.762 1.281 1.270];
        k_exp_IAEA_CRP = [0.04691 0.04575 0.04785 0.04579 0.03850 0.03817];
    elseif strcmp(temperature, '900C')
        time_exp_IAEA_CRP = [30 30 120 120 480 480]/86400;
        wg_exp_IAEA_CRP = [0.579 0.556 0.998 1.017 0.503 0.448];
        k_exp_IAEA_CRP = [0.10572 0.10148 0.09107 0.09286 0.02296 0.02046];
    end
end


%% Thomazet
if strcmp(alloy, 'Zr10Nb') && strcmp(temperature, '360C')
    time_fh_Thomazet_M5 = [0 1 2 25 50 75 100 150 200 250 325 345 355 375 390 410 430 455 525 550 575 590 615 640];
    fh_Thomazet_M5 = [1.35 2.1 2.1 2.1 2.3 2.5 2.05 2.35 2.75 3.55 3.85 4.4 4.9 5.25 5.3 5.7 5.75 5.8 5.85 5.9 6.3 6.7 6.6 7.5];
    time_fh_Thomazet_Zry4 = [0 1 2.5 4 24.7 50 70 75 80 85 90 101 125 152 162 167 172 177 182 198];
    fh_Thomazet_Zry4 = [1 1.4 1.63 1.7 2.05 2.1 2.52 2.86 2.95 2.97 3.13 3.39 3.635 4 3.95 4.25 4 4.27 4.41 4.9];
 
    time_d_exp_Thomazet_M5 = [0 25 50 75 100 135 150 200 215 250 300 325 350 370 380 395 420 445 475 500 520 530 540 555 570 585 600 615 640 660 690 715 740 775 805 835 860 890 925 950 980 1015 1050 1080 1110 1140 1175 1200 1225 1255 1285 1315 1330 1360 1400 1425 1450 1485 1515 1550 1575];
    d_exp_Thomazet_M5 = [0 0.95 1.45 1.65 2 2.5 2.75 3 3.25 3.5 4.2 4.25 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7 7.25 7.3 7.5 7.75 8 8.25 8.5 8.75 9.25 9.25 9.5 10 10.25 10.5 10.75 11.25 11.75 12 12.5 12.75 13 13.25 13.5 13.75 14.25 14.5 15 15.5 15.75 16 16.25 16.75 17 17.25 17.25 17.75 18.25 18.5 18.75];
    time_d_exp_Thomazet_Zry4 = [1 3 5 10 15 20 25 39 50 71 76 81 86 90 96 101 106 111 126 131 136 152 162 167 172 177 182 200 215 250 300 325 350 370 380 400 425 450 475 500 525 530 540 555 570 585 600 615 630 655 690 710 740 770 800 830 855 885 915 950 975 1010 1050 1075 1110 1140 1165 1200 1225 1250 1280 1315 1325 1360 1390 1420 1435 1465 1500 1530 1560];
    d_exp_Thomazet_Zry4 = [0.45 0.6 0.7 1 1.1 1.2 1.25 1.45 1.6 1.85 2.1 2.25 2.35 2.4 2.55 2.7 2.9 3.05 3.25 3.3 3.35 3.45 3.5 3.7 3.95 4.2 4.3 5 5.1 5.6 7 7.5 8.25 8.75 9 9.25 9.75 10.5 11 12 12.75 13 13 13.5 13.75 14.5 14.5 15 15.25 16 16.5 17 18 18.5 19 20 20.5 21.5 22.25 22.75 23.5 24.25 25.25 26 26.75 27.5 28 28.75 29.5 30.25 30.75 31.5 31.75 32.5 33.5 34 34 35 36 36.5 37];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% High Temperature experimental data

%% Experiment UW-Madison

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '950C')
        time_exp_UW = [500 1000 1500];       % [s]
        d_exp_UW = [];
        alpha_exp_UW = [];
        wg_exp_UW = [3.92 5.57 6.01];     % [mg/cm^2]
    end
    if strcmp(temperature, '1100C')
        time_exp_UW = [100 200 500 1000];% 500 1300 1500];       % [s]
        d_exp_UW = [22 30 47 64];
        alpha_exp_UW = [20 31 48 63];
        wg_exp_UW = [3.99 5.83 9.33 12.50];% 9.52 14.80 15.92];     % [mg/cm^2]
    end
    if strcmp(temperature, '1200C')
        time_exp_UW = [50 200 500 1000 1500];% 500 800 1300];               % [s]
        d_exp_UW = [24 45 68 94 117];
        alpha_exp_UW = [28 59 88 123 147];
        wg_exp_UW = [4.90 10.19 15.45 21.34 25.91];% 14.85 15.81 22.52];       % [mg/cm^2]
    end
end


%% Experiment multiple oxidations MIT

if strcmp(alloy, 'Zry4') && strcmp(temperature, '1200C')
    time_exp_MIT = [300 600 900];
    d_exp_MIT = [35.9 52 71.3];
    alpha_exp_MIT = [41.3 58.3 81];
    wg_exp_MIT = [NaN 16.3 19.8];
end


%% Cathcart-Pawel experiment

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '950C')
        time_exp_CP = [1344.9 977.9 1595.7 690.6 774.6 1300.2 948.7 1634.5 412.2 428.7];% [s]
        d_exp_CP = [26.5 23.1 26.2 21.3 22.5 25.1 22.0 26.3 17.9 18.0];       % [um]
        alpha_exp_CP = [24.1 18.9 25.2 15.7 16.9 22.3 19.9 26.4 13.0 15.3];   % [um]
        wg_exp_CP = [4.62 3.97 4.61 3.62 3.83 4.37 3.83 4.65 3.04 3.11];      % [mg/cm^2]
        time_CP_max = 1500;
    end
    if strcmp(temperature, '1000C')
        time_exp_CP = [1130.5 725.2 301.1 925.0 300.2 483.4 678.6 1104.1 915.5 475.7 493.4 299.4];      % [s]
        d_exp_CP = [40.5 32.6 22.5 37.4 20.9 26.5 29.9 38.6 34.5 25.9 27.1 21.3];       % [um]
        alpha_exp_CP = [23.7 23.1 15.2 24.6 14.4 18.0 20.8 26.1 24.6 17.6 19.5 13.8];   % [um]
        wg_exp_CP = [7.04 5.77 3.95 6.57 3.69 4.67 5.30 6.81 6.13 4.57 4.80 3.73];      % [mg/cm^2]
        time_exp_CP_MAXI = [367.4 866.5 99.4 640.5 234.4 512.4];    %[s]
        d_exp_CP_MAXI = [25.8 40.6 13.7 35.3 20.6 31.0];            % [um]
        alpha_exp_CP_MAXI = [17.6 26.0 10.2 22.6 14.1 19.3];        % [um]
        wg_exp_CP_MAXI = [4.547 7.104 2.436 6.173 3.632 5.410];     % [um]
        time_CP_max = 1000;
    end
    if strcmp(temperature, '1100C')
        time_exp_CP = [495.6 252.7 313.4 512.4 348.1 283.4 252.1 159.5 173.9 449.3 458.9];
        d_exp_CP = [45.7 33.5 36.6 47.1 37.2 40.8 33.3 26.0 29.4 45.7 44.6];
        alpha_exp_CP = [41.6 33.0 36.6 45.4 37.0 37.0 24.5 23.1 26.9 43.5 44.5];
        wg_exp_CP = [8.77 6.49 7.12 9.10 7.25 7.82 6.22 4.97 5.61 8.79 8.66];
        time_CP_max = 500;
    end
    if strcmp(temperature, '1200C')
        time_exp_CP = [236.4 160.2 126.6 280.2 234.6 56.1 171.8 66.1 111.2];
        d_exp_CP = [49.4 42.8 38.9 53.5 50.2 26.6 43.2 28.6 37.2];
        alpha_exp_CP = [53.5 44.7 41.0 54.9 51.6 28.2 44.0 30.9 39.7];
        wg_exp_CP = [10.03 8.59 7.81 10.78 10.09 5.33 8.66 5.75 7.46];
        time_CP_max = 300;
    end
    if strcmp(temperature, '1300C')
        time_exp_CP = [138.4 151.4 83.4 124.4 123.2 32.9 30.0 111.1 62.2 79.9 59.3 57.7 34.3];
        d_exp_CP = [56.2 59.7 44.9 51.7 53.1 28.9 27.9 50.1 38.2 41.2 36.4 36.5 28.0];
        alpha_exp_CP = [66.9 71.8 53.4 61.7 60.6 32.5 31.7 58.9 45.2 48.6 42.8 42.2 32.4];
        wg_exp_CP = [11.87 12.61 9.45 10.96 11.13 6.01 5.81 10.57 8.05 8.72 7.68 7.67 5.89];
        time_CP_max = 200;
    end
    if strcmp(temperature, '1400C')
        time_exp_CP = [72.7 62.2 29.7 28.6 45.7 42.6 14.3 13.5 11.2 10.4 56.4 54.2 44.8 43.4 27.5 29.8 74.1 78.9 29.3 26.5];
        d_exp_CP = [58.4 53.4 38.7 35.2 45.7 44.7 28.0 26.1 24.5 23.4 51.6 48.8 47.1 45.1 39.1 39.2 59.0 60.3 39.0 35.1];
        alpha_exp_CP = [78.3 70.9 51.0 49.2 62.9 60.4 36.8 35.7 31.8 30.7 67.5 65.7 62.2 60.4 51.1 51.4 77.3 77.6 50.4 47.2];
        wg_exp_CP = [12.84 11.74 8.45 7.84 10.12 9.85 6.08 5.74 5.31 5.09 11.29 10.79 10.29 9.92 8.47 8.53 12.92 13.18 8.46 7.73];
        time_CP_max = 100;
    end
    if strcmp(temperature, '1500C')
        time_exp_CP = [31.8 28.9 33.0 31.2 47.0 42.6 22.8 21.9 7.6 8.2 53.2 9.2 8.9 49.0 50.0 42.4 37.8 13.9 13.2];
        d_exp_CP = [53.6 50.7 57.8 50.3 65.6 61.8 45.8 43.8 27.3 28.5 73.7 30.0 28.9 62.6 62.5 65.9 57.5 37.5 35.5];
        alpha_exp_CP = [72.1 68.1 75.7 74.2 88.2 83.8 63.6 62.8 38.7 38.5 97.6 40.7 40.8 93.9 91.4 85.6 80.3 48.3 48.6];
        wg_exp_CP = [11.98 11.34 12.76 11.53 14.65 13.84 10.29 9.93 6.13 6.34 16.29 6.68 6.50 14.42 14.34 14.51 12.97 8.25 7.93];
        time_CP_max = 50;
    end
end


%% Sawarn experiment

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '1000C')
        time_exp_Sawarn = [60 120 300 600 900];                 % [s]
        d_exp_Sawarn = [25.11 37.08 39.59 50.46 67.49];         % [um]
        d_exp_Sawarn_error = [1.59 2.67 1.87 3.47 2.10];        % [um]
        alpha_exp_Sawarn = [25.4 41.71 70.52 91.83 111.06];         % [um]
        alpha_exp_Sawarn_error = [6.51 14.34 24.00 32.48 26.99];% [um]
    end
    if strcmp(temperature, '1100C')
        time_exp_Sawarn = [60 120 300 600 900];                 % [s]
        d_exp_Sawarn = [49.6 56.35 75.7 109.85 133.74];         % [um]
        d_exp_Sawarn_error = [1.86 4.91 7.65 2.87 5.47];        % [um]
        alpha_exp_Sawarn = [36.76 52 83.3 139.11 163.55];           % [um]
        alpha_exp_Sawarn_error = [1.84 4.36 21.49 8.84 7.31];   % [um]
    end
    if strcmp(temperature, '1200C')
        time_exp_Sawarn = [60 120 300 600 900];                 % [s]
        d_exp_Sawarn = [62.26 82.8 126.94 179.89 198.98];       % [um]
        d_exp_Sawarn_error = [4.05 2.37 17.49 4.86 6.45];       % [um]
        alpha_exp_Sawarn = [77 116.8 180.2 305.1 311.26];           % [um]
        alpha_exp_Sawarn_error = [14.10 6.68 43.73 9.43 9.34];   % [um]
    end
end


%% NRC experiment

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '1000C')
        time_exp_NRC = 3364;    % [s]
        d_exp_NRC = 82.5;
        wg_exp_NRC = 14.6;      % [mg/cm^2]
    end
    if strcmp(temperature, '1100C')
        time_exp_NRC = 1065;
        d_exp_NRC = 69;
        wg_exp_NRC = 13.2;      % [mg/cm^2]
    end
    if strcmp(temperature, '1200C')
        time_exp_NRC = [166 400];
        d_exp_NRC = [41.5 67];
        alpha_exp_NRC = [48.7 108.8];
        wg_exp_NRC = [8.35 13.5];       % [mg/cm^2]
    end
end


%% Lemmon experiment

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '1000C')
        time_exp_Lemmon = [120 7200];
        d_exp_Lemmon = [24 230];
        alpha_exp_Lemmon = [8 200];
        wg_exp_Lemmon = Mo / VV * [5.2 31];
    end
    if strcmp(temperature, '1200C')
        time_exp_Lemmon = [120 3600];
        d_exp_Lemmon = [30 220];
        alpha_exp_Lemmon = [28 45];
        wg_exp_Lemmon = Mo / VV * [8.8 74];
    end
    if strcmp(temperature, '1400C')
        time_exp_Lemmon = [120 2100];
        d_exp_Lemmon = [60 280];
        alpha_exp_Lemmon = [84 380];
        wg_exp_Lemmon = Mo / VV * [20 82];
    end
end


%% Biederman (WPI) experiment

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '1000C')
        time_exp_WPI = [100 100 100 500 500 500];
        wg_exp_WPI = [2.8 2.2 1.95 4.0 4.3 3.66];   % [g/cm^2]
    end
    if strcmp(temperature, '1100C')
        time_exp_WPI = [100 100 100 500 500 500];
        wg_exp_WPI = [5.7 4.6 4.19 7.6 8.4 7.76];   % [g/cm^2]
    end
    if strcmp(temperature, '1200C')
        time_exp_WPI = [100 100 100 500 500 500];
        wg_exp_WPI = [8.1 6.4 6.41 13.9 13.1 12.3];   % [g/cm^2]
    end
    if strcmp(temperature, '1250C')
        time_exp_WPI = [50 50 50 100 100 100];
        wg_exp_WPI = [5.9 5.7 5.02 8.9 7.0 6.70];   % [g/cm^2]
    end
end

%% Brachet experimental data

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '1200C')
        time_wg_exp_Brachet_Zry4 = [54.8 196.0 506.3 529.0 1482.3 1513.2 2116 3249 3249];
        wg_exp_Brachet_Zry4 = [4.1 8.2 13.5  14.3 22.3 23.3 26 31.7 33.3];   % [um]
        time_wg_exp_Brachet_M5 = [56.3 196 449.4 506.3 552.3 552.3 1444 1444 3025 3025];
        wg_exp_Brachet_M5 = [4 8.5 14.1 13.5 14.9 16 23 23.8 35 36];   % [um]
        
        time_d_exp_Brachet_Zry4 = [56.3 190.4 519.8 1497.7];
        d_exp_Brachet_Zry4 = [23 43 71 113];   % [um]
        time_d_exp_Brachet_M5 = [60.8 201.6 566.4 1421.3];
        d_exp_Brachet_M5 = [18.5 45 73 122];   % [um]
        
        time_alpha_exp_Brachet_Zry4 = [56.3 193.2 524.5 1505.4];
        alpha_exp_Brachet_Zry4 = [20 45.5 86 140];   % [um]
        time_alpha_exp_Brachet_M5 = [59.3 198.8 557.0 1413.8];
        alpha_exp_Brachet_M5 = [20 51 82 148];   % [um]
    end
end


%% ORNL experiment (Ben Garrison)

if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '1000C')
        time_exp_ORNL = [990 1890 2850 4000];
        wg_exp_ORNL = [6.05 8.10 8.60 9.93];            % [mg/cm^2]
        d_exp_ORNL = [36.30 43.88 50.57 58.67];         % [um]
    end
    if strcmp(temperature, '1100C')
        time_exp_ORNL = [200 380 700 950];
        wg_exp_ORNL = [6.38 8.16 10.88 12.19];          % [mg/cm^2]
        d_exp_ORNL = [36.77 42.84 57.17 64.01];         % [um]
    end
    if strcmp(temperature, '1200C')
        time_exp_ORNL = [46 126 271 406];
        wg_exp_ORNL = [7.48 9.17 12.35 14.48];          % [mg/cm^2]
        d_exp_ORNL = [32.50 41.70 54.17 64.30];         % [um]
    end
end


%% OAH-ABA
if strcmp(alloy, 'Zry4')
    if strcmp(temperature, '900C')
        time_exp_OAH_ABA_Zr10Nb = [360 1000 3000 7000 11000 14000]; % [s]
        wg_exp_OAH_ABA_Zr10Nb = [1.208 2.767 6.455 10.297 14.178 14.550]; % [mg/cm^2]
        k_exp_OAH_ABA_Zr10Nb = [0.064 0.087 0.118 0.123 0.135 0.123]; % [mg/cm^2*s^0.5]
        
        time_exp_OAH_ABA_Zry4 = [300 1000 5000 11000]; % [s]
        wg_exp_OAH_ABA_Zry4 = [1.72 2.64 4.50 5.32]; % [mg/cm^2]
        k_exp_OAH_ABA_Zry4 = [0.100 0.084 0.064 0.051]; % [mg/cm^2*s^0.5]
    elseif strcmp(temperature, '1000C')
        time_exp_OAH_ABA_Zr10Nb = [100 700 1200 1800 3600 6000]; % [s]
        wg_exp_OAH_ABA_Zr10Nb = [1.454 4.637 7.100 12.710 17.764 22.613]; % [mg/cm^2]
        k_exp_OAH_ABA_Zr10Nb = [0.145 0.175 0.205 0.300 0.296 0.292]; % [mg/cm^2*s^0.5]
        
        time_exp_OAH_ABA_Zry4 = [87 464 2600 3300 4090 7270 11360]; % [s]
        wg_exp_OAH_ABA_Zry4 = [2.19 4.56 9.25 11.52 15.42 33.40 58.61]; % [mg/cm^2]
        k_exp_OAH_ABA_Zry4 = [0.235 0.212 0.181 0.201 0.241 0.392 0.550]; % [mg/cm^2*s^0.5]
    elseif strcmp(temperature, '1100C')
        time_exp_OAH_ABA_Zr10Nb = [19 133 704 1500 2400 5000]; % [s]
        wg_exp_OAH_ABA_Zr10Nb = [1.166 3.533 8.548 12.555 16.305 23.554]; % [mg/cm^2]
        k_exp_OAH_ABA_Zr10Nb = [0.268 0.306 0.322 0.324 0.333 0.333]; % [mg/cm^2*s^0.5]
        
        time_exp_OAH_ABA_Zry4 = [27 102 398 900 1500 3000]; % [s]
        wg_exp_OAH_ABA_Zry4 = [2.14 4.14 7.64 11.54 14.78 20.25]; % [mg/cm^2]
        k_exp_OAH_ABA_Zry4 = [0.412 0.410 0.383 0.385 0.382 0.370]; % [mg/cm^2*s^0.5]
    elseif strcmp(temperature, '1200C')
        time_exp_OAH_ABA_Zr10Nb = [7 49 167 380 646 1205]; % [s]
        wg_exp_OAH_ABA_Zr10Nb = [1.641 3.693 7.620 11.035 14.455 19.783]; % [mg/cm^2]
        k_exp_OAH_ABA_Zr10Nb = [0.620 0.528 0.590 0.566 0.569 0.570]; % [mg/cm^2*s^0.5]
        
        time_exp_OAH_ABA_Zry4 = [10 40 163 367 790 1100]; % [s]
        wg_exp_OAH_ABA_Zry4 = [2.63 4.35 7.93 11.63 16.54 19.29]; % [mg/cm^2]
        k_exp_OAH_ABA_Zry4 = [0.832 0.688 0.621 0.607 0.588 0.582]; % [mg/cm^2*s^0.5]
    end
end


%% Empirical models - time range isotherm

if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
        time_CP = 0:dt:time_CP_max;
        time_CP_extended = 0:dt:1500;
        time_BJ = 0:dt:3000;
        time_Lemmon = 0:dt:3000;
        time_WPI = 0:dt:1000;
        time_Kawasaki = 0:dt:1500;
        time_Leistikow = 0:dt:1800;
        time_Urbanic = 0:dt:2000;
    end
end


%% Empirical models - weight gain isotherm

if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
        wg_CP = sqrt(2 * 0.1811 * exp(-39940/(RR*T(1))) .* time_CP) * 1e3;     % [mg/cm^2]
        wg_CP_extended = sqrt(2 * 0.1811 * exp(-39940/(RR*T(1))) .* time_CP_extended) * 1e3;     % [mg/cm^2]
%         wg_CP_grad = sqrt(2 * 0.1680 * exp(-39870/(RR*T(1))) .* time_CP) * 1e3;    % [mg/cm^2]
        wg_BJ = sqrt((2 * Mo / MZr)^2 * 33.3e6 * exp(-45500/(RR*T(1))) .* time_BJ);     % [mg/cm^2]
%         wg_Baker = sqrt(2 * 2.0496 * exp(-45500/(RR*T(1))) .* time_BJ) * 1e3;     % [mg/cm^2]
        wg_Lemmon = sqrt((Mo / 22.71)^2 * 0.1132e6 * exp(-34000/(RR*T(1))) .* time_Lemmon);    % [mg/cm^2]
%         wg_Lemm = sqrt(2 * 0.028875 * exp(-34000/(RR*T(1))) .* time_Lemmon) * 1e3;    % [mg/cm^2]
        wg_Leistikow = sqrt(2 * 2.142 * exp(-47640/(RR*T(1))) .* time_Leistikow) * 1e3;  % [mg/cm^2]
        wg_WPI = 195.3 * sqrt(exp(-33370/(RR*T(1))) .* time_WPI);            % [mg/cm^2]
%         wg_Biederman = sqrt(2 * 0.01907 * exp(-33370/(RR*T(1))) .* time_WPI) * 1e3;      % [mg/cm^2]
        wg_Kawasaki = sqrt(0.468 * exp(-40710/(RR*T(1))) .* time_Kawasaki) * 1e3;   % [mg/cm^2]
%         wg_Kawa = sqrt(2 * 0.1994 * exp(-40500/(RR*T(1))) .* time_Kawasaki) * 1e3;  % [mg/cm^2]
        wg_Urbanic = sqrt((2 * Mo / MZr)^2 * 2.96e5 * exp(-33421/(RR*T(1))) .* time_Urbanic);     % [mg/cm^2] (below 1580C)
        
        wg_Hobson = sqrt(2 * 0.1553 * exp(-39290/(RR*T(1))) .* time) * 1e3;    % [mg/cm^2]
%         wg_Hobson_point = sqrt(2 * 0.008311 * exp(-31110/(RR*T(1))) .* time) * 1e3;    % [mg/cm^2]
        wg_Klepfer = sqrt(2 * 0.02203 * exp(-33500/(RR*T(1))) .* time) * 1e3;  % [mg/cm^2]
        wg_Parsons = sqrt( (2 * Mo / MZr)^2 * 25.14e6 * exp(-45600/(RR*T(1))) .* time); % [mg/cm^2]
    end
end


%% Empirical models - oxide thickness isotherm

if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
        d_CP = sqrt(2 * 0.01126 * exp(-35890/(RR*T(1))) .* time_CP);                      % [cm]
        d_CP_extended = sqrt(2 * 0.01126 * exp(-35890/(RR*T(1))) .* time_CP_extended);    % [cm]
        d_BJ = wg_BJ ./ (14.77 * 1e2);                                                    % [cm]
        d_WPI = 364.4 * sqrt(exp(-27840/(RR*T(1))) .* time_WPI) * 1e-4;                   % [cm]
        d_Leistikow = sqrt(7.82e-2 * exp(-168100/(R*T(1))) .* time_Leistikow);            % [cm]
        d_Kawasaki = sqrt(0.0215 * exp(-35860/(RR*T(1))) .* time_Kawasaki);               % [cm]
        d_Urbanic = sqrt(0.036^2 * exp(-26995/(RR*T(1))) .* time_Urbanic);                % [cm]
    end
end


%% Empirical models - alpha thickness isotherm

if strcmp(mode, 'Isotherm')
    if strcmp(alloy, 'Zry4') && T0 > T_transf
        alpha_CP = sqrt(2 * 0.7615 * exp(-48140/(RR*T(1))) .* time_CP);                         % [cm]
        alpha_CP_extended = sqrt(2 * 0.7615 * exp(-48140/(RR*T(1))) .* time_CP_extended);       % [cm]
        alpha_WPI = 2257.5 * sqrt(exp(-37690/(RR*T(1))) .* time_WPI) * 1e-4;                    % [cm]
        alpha_Kawasaki = sqrt(0.306 * exp(-38970/(RR*T(1))) .* time_Kawasaki) - d_Kawasaki;     % [cm]
        alpha_Urbanic = sqrt(0.390^2 * exp(-39402/(RR*T(1))) .* time_Urbanic);                  % [cm]
    end
end


%% Import values from BISON simulation
%{
if strcmp(alloy, 'Zry4') && strcmp(temperature, '1200C')
    Cathcart_file = "../../../BISON/inputs/oxidation/bison_Cathcart_1200C_500s.csv";
    Leistikow_file = "../../../BISON/inputs/oxidation/bison_Leistikow_1200C_500s.csv";
    C4_isotherm_file = "../../../BISON/inputs/oxidation/C4_isotherm.csv";
    C4_O_file = "../../../BISON/inputs/oxidation/C4_O_1200C_500s.csv";
%     C4_O_file_fine = "../../../BISON/inputs/oxidation/C4_O_1200C_500s_Tcste.csv";
    
    data_bison_Cathcart = csvread(Cathcart_file, 1, 0);
    data_bison_Leistikow = csvread(Leistikow_file, 1, 0);
    data_bison_isotherm = csvread(C4_isotherm_file, 1, 0);
    data_bison_C4_O = csvread(C4_O_file, 1, 0);
%     data_bison_C4_O_fine = csvread(C4_O_file_fine, 1, 0);
    
    t_bison_Cathcart = data_bison_Cathcart(:, 59);      % [s]
    d_bison_Cathcart = data_bison_Cathcart(:, 14);      % [m]
    wg_bison_Cathcart = data_bison_Cathcart(:, 21);     % [kg/m^2]
    t_bison_Leistikow = data_bison_Leistikow(:, 59);    % [s]
    d_bison_Leistikow = data_bison_Leistikow(:, 14);    % [m]
    wg_bison_Leistikow = data_bison_Leistikow(:, 21);   % [kg/m^2]
    t_bison_isotherm = data_bison_isotherm(:, 59);      % [s]
    d_bison_isotherm = data_bison_isotherm(:, 14);      % [m]
    wg_bison_isotherm = data_bison_isotherm(:, 21);     % [kg/m^2]
    t_bison_C4_O = data_bison_C4_O(:, 59);              % [s]
    d_bison_C4_O = data_bison_C4_O(:, 14);              % [m]
    wg_bison_C4_O = data_bison_C4_O(:, 21);             % [kg/m^2]
%     t_bison_C4_O_fine = data_bison_C4_O_fine(:, 59);              % [s]
%     d_bison_C4_O_fine = data_bison_C4_O_fine(:, 14);              % [m]
%     wg_bison_C4_O_fine = data_bison_C4_O_fine(:, 21);             % [kg/m^2]
end
%}


%% Fitting of C4 model to parabolic law

fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf('R-square values:\n\n');

[fit_wg, fit_wg_stat] = fit(time(2:N+1)', wg(2:N+1)', 'power1', 'Lower', [-inf 0.5], 'Upper', [inf 0.5]);
fit_wg_coef = coeffvalues(fit_wg);
wg_fit = fit_wg_coef(1) * time.^(0.5);

fprintf('weight gain R-square: %.5f\n', fit_wg_stat.rsquare);

[fit_d, fit_d_stat] = fit(time(2:N+1)', d_protect(2:N+1)', 'power1', 'Lower', [-inf 0.5], 'Upper', [inf 0.5]);
fit_d_coef = coeffvalues(fit_d);
d_fit = fit_d_coef(1) * time.^(0.5);

fprintf('oxide thickness R-square: %.5f\n', fit_d_stat.rsquare);

[fit_alpha, fit_alpha_stat] = fit(time(2:N+1)', alpha(2:N+1)', 'power1', 'Lower', [-inf 0.5], 'Upper', [inf 0.5]);
fit_alpha_coef = coeffvalues(fit_alpha);
alpha_fit = fit_alpha_coef(1) * time.^(0.5);

fprintf('alpha thickness R-square: %.5f\n', fit_alpha_stat.rsquare);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Temperature transient experimental data

%% Empirical models - weight gain transient

if strcmp(mode, 'Linear') || strcmp(mode, 'LOCA')
    wg_CP = zeros(1, N+1);
    wg_BJ = zeros(1, N+1);
    
    for n = 2:N+1
        wg_CP(n) = sqrt((wg_CP(n-1) * 1e-3)^2 + 2 * 0.1811 * exp(-39940/(RR*T(n))) * dt) * 1e3;
        wg_BJ(n) = sqrt(wg_BJ(n-1)^2 + (2 * Mo / MZr)^2 * 33.3e6 * exp(-45500/(RR*T(n))) * dt);
    end
end


%% Empirical models - oxide thickness transient

if strcmp(mode, 'Linear') || strcmp(mode, 'LOCA')
    d_CP = zeros(1, N+1);
    d_BJ = wg_BJ ./ (14.77 * 1e2);
    
    for n = 2:N+1
        d_CP(n) = sqrt(d_CP(n-1)^2 + 2 * 0.01126 * exp(-35890/(RR*T(n))) * dt);
    end
end

%% Empirical models - alpha thickness transient

if strcmp(mode, 'Linear') || strcmp(mode, 'LOCA')
    alpha_CP = zeros(1, N+1);
    
    for n = 2:N+1
        alpha_CP(n) = sqrt(alpha_CP(n-1)^2 + 2 * 0.7615 * exp(-48140/(RR*T(n))) * dt);
    end
end


%% Linear CINOG fit (700-1200C)

if strcmp(alloy, 'Zry4') && strcmp(temperature, '700-1200C') && strcmp(exposure_time, '500s')
    wg_CINOG = 4.75e-8 * time.^3 - 6.97e-6 * time.^2 + 2.66e-3 * time + 1.08e-1;
    wg_Couet = 1.74e-8 * time.^3 - 2.4e-7 * time.^2 + 2.45e-3 * time;
    wg_Couetpp = 2.48e-8 * time.^3 - 3.32e-6 * time.^2 + 3.27e-3 * time;
    wg_BJ_CINOG = 6.56e-8 * time.^3 - 7.42e-6 * time.^2 + 2.69e-3 * time + 8.52e-2;
elseif strcmp(alloy, 'Zry4') && strcmp(temperature, '700-1200C') && strcmp(exposure_time, '5000s')
    wg_CINOG = 1.5e-10 * time.^3 - 2.55e-7 * time.^2 + 1.10e-3 * time - 1.61e-2;
    wg_Couet = 5.96e-11 * time.^3 - 3.73e-8 * time.^2 + 8.26e-4 * time;
    wg_Couetpp = 7.91e-11 * time.^3 - 8.37e-8 * time.^2 + 9.25e-4 * time;
    wg_BJ_CINOG = 2.19e-10 * time.^3 - 3.33e-7 * time.^2 + 1.15e-3 * time - 1.54e-2;
end

%% CINOG (Zanella report)

if strcmp(alloy, 'Zry4') && strcmp(temperature, 'CINOG1')
    wg_exp_CINOG = 5.3;
    wg_CINOG = 7.0;
elseif strcmp(alloy, 'Zry4') && strcmp(temperature, 'CINOG2')
    wg_exp_CINOG = 2.5;
    wg_CINOG = 2.2;
elseif strcmp(alloy, 'Zry4') && strcmp(temperature, 'CINOG3')
    wg_exp_CINOG = 18.2;
    wg_CINOG = 18.0;
elseif strcmp(alloy, 'Zry4') && strcmp(temperature, 'CINOG4')
    wg_exp_CINOG = 6.1;
    wg_CINOG = 5.7;
elseif strcmp(alloy, 'Zry4') && strcmp(temperature, 'CINOG5')
    wg_exp_CINOG = 16.8;
    wg_CINOG = 16.6;
elseif strcmp(alloy, 'Zry4') && strcmp(temperature, 'CINOG6')
    wg_exp_CINOG = 8.2;
    wg_CINOG = 5.3;
end


%% ANL (LOCA)
if strcmp(mode, 'LOCA')
    if strcmp(alloy, 'Zry4') % Exposure Times here are the full experiment time
        if strcmp(temperature, '1200C') && strcmp(exposure_time, '600s')
            wg_CP_ANL = 12.3;
            wg_exp_ANL = 14.9;
        elseif  strcmp(temperature, '1200C') && strcmp(exposure_time, '900s')
            wg_CP_ANL = 16.6;
            wg_exp_ANL = 19.9;
        elseif  strcmp(temperature, '1200C') && strcmp(exposure_time, '1500s')
            wg_CP_ANL = 23.7;
            wg_exp_ANL = 26.1;
        elseif  strcmp(temperature, '1100C') && strcmp(exposure_time, '2100s')
            wg_CP_ANL = 17.0;
            wg_exp_ANL = 20.4;
        end
    end
end


%% Leistikow Peak

if strcmp(mode, 'LOCA')
    if strcmp(alloy, 'Zry4')
        if strcmp(temperature, 'Peak950C')
            wg_exp_Leistikow = [61.9 49.3 46.2] * 1e-2;         % [mg/cm^2]
        elseif  strcmp(temperature, 'Peak1000C')
            wg_exp_Leistikow = [56.4 73.1 65.4 47.3] * 1e-2;    % [mg/cm^2]
        elseif  strcmp(temperature, 'Peak1100C')
            wg_exp_Leistikow = [119 122 115] * 1e-2;            % [mg/cm^2]
        elseif  strcmp(temperature, 'Peak1200C')
            wg_exp_Leistikow = [203 191 191] * 1e-2;            % [mg/cm^2]
        end
    end
end


%% Leistikow LOCA
if strcmp(mode, 'LOCA')
    if strcmp(alloy, 'Zry4')
        if strcmp(temperature, 'LeistikowLOCA')
            time_Leistikow_LOCA = 0:5:200;
            wg_Leistikow_LOCA = [0 0.1 0.3 0.38 0.42 0.42 0.44 0.46 0.48 0.5 0.55 0.62 0.69 0.75 0.86 1 1.15 1.3 1.47 1.74 2 2.25 2.51 2.8 3.08 3.35 3.6 3.82 4.04 4.25 4.44 4.62 4.8 4.98 5.15 5.3 5.47 5.51 5.52 5.52 5.52];
            d_Leistikow_LOCA = [0 0.5 1.6 2.1 2.4 2.4 2.4 2.6 2.8 2.9 3.2 3.6 4 4.3 4.9 5.6 6.2 7 7.9 9 10.1 11.3 12.5 13.8 15 16.2 17.4 18.4 19.4 20.4 21.3 22.2 23 24 24.7 25.5 26.3 26.3 26.2 26.1 26];
            alpha_Leistikow_LOCA = [0 0.5 1.6 2.1 2.4 2.4 2.4 2.4 2.4 2.4 2.7 3.1 3.5 3.8 4.5 5.6 6.5 7.5 9 11 13.1 15.3 17.4 20 22.5 24.9 27 28.9 30.8 32.5 34 35.5 37 38.4 39.7 40.9 42 53 55.6 55.6 55.6]; %microns
            
            time_exp_Leistikow_LOCA = [24 24 72 72 91 92 110 110 203 203];
            wg_exp_Leistikow_LOCA = [0.48 0.61 0.74 0.82 1.44 1.48 2.3 2.88 6.2 5.4];
            d_exp_Leistikow_LOCA = 44;
            alpha_exp_Leistikow_LOCA = 34;
        end
    end
end
