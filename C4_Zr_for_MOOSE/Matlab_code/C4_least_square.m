%% C4_least_square
%%Determine a least square criteria to compare the models and experimental results
%%The criteria is the normalized root-mean-square deviation

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

% Create an array that regroup all the data

if strcmp(alloy, 'Zry4') && T0 > T_transf
    N_model = 7;
    N_time_array = max([length(time), length(time_CP), length(time_BJ), length(time_Leistikow), length(time_WPI), length(time_Kawasaki), length(time_Urbanic)]);
    
    time_model_array = nan(N_time_array, N_model);
    time_model_array(1:length(time),1) = time';
    time_model_array(1:length(time_CP),2) = time_CP';
    time_model_array(1:length(time_BJ),3) = time_BJ';
    time_model_array(1:length(time_Leistikow),4) = time_Leistikow';
    time_model_array(1:length(time_WPI),5) = time_WPI';
    time_model_array(1:length(time_Kawasaki),6) = time_Kawasaki';
    time_model_array(1:length(time_Urbanic),7) = time_Urbanic';
    
    % Model order: C4 ; CP ; BJ ; Leistikow ; WPI ; Kawasaki
    wg_model_array = nan(N_time_array, N_model);
    wg_model_array(1:length(time),1) = wg';
    wg_model_array(1:length(time_CP),2) = wg_CP';
    wg_model_array(1:length(time_BJ),3) = wg_BJ';
    wg_model_array(1:length(time_Leistikow),4) = wg_Leistikow';
    wg_model_array(1:length(time_WPI),5) = wg_WPI';
    wg_model_array(1:length(time_Kawasaki),6) = wg_Kawasaki';
    wg_model_array(1:length(time_Urbanic),7) = wg_Urbanic';
    
    % Model order: C4 ; CP ; BJ ; Leistikow ; WPI ; Kawasaki ; Urbanic
    d_model_array = nan(N_time_array, N_model);
    d_model_array(1:length(time),1) = d_protect';
    d_model_array(1:length(time_CP),2) = d_CP';
    d_model_array(1:length(time_BJ),3) = d_BJ';
    d_model_array(1:length(time_Leistikow),4) = d_Leistikow';
    d_model_array(1:length(time_WPI),5) = d_WPI';
    d_model_array(1:length(time_Kawasaki),6) = d_Kawasaki';
    d_model_array(1:length(time_Urbanic),7) = d_Urbanic';
    
    % Model order: C4 ; CP ; WPI ; Kawasaki ; Urbanic
    alpha_model_array = nan(N_time_array, N_model);
    alpha_model_array(1:length(time),1) = alpha';
    alpha_model_array(1:length(time_CP),2) = alpha_CP';
    alpha_model_array(1:length(time_WPI),5) = alpha_WPI';
    alpha_model_array(1:length(time_Kawasaki),6) = alpha_Kawasaki';
    alpha_model_array(1:length(time_Urbanic),7) = alpha_Urbanic';
    
    
    wg_diff = zeros(N_model);
    d_diff = zeros(N_model);
    alpha_diff = zeros(N_model);
    
    for i = 1:N_model
        for j = 1:N_model
            if i == j
                wg_diff(i,j) = 1;
                d_diff(i,j) = 1;
                alpha_diff(i,j) = 1;
            else
                for k = 3:N_time_array
                    if (isnan(wg_model_array(k,i)) || isnan(wg_model_array(k,j)))
                    else
                        wg_diff(i,j) = sqrt(wg_diff(i,j)^2 + ((wg_model_array(k,i) - wg_model_array(k,j)) / wg_model_array(k,j))^2 / N);
                    end
                    
                    if (isnan(d_model_array(k,i)) || isnan(d_model_array(k,j)))
                    else
                        d_diff(i,j) = sqrt(d_diff(i,j)^2 + ((d_model_array(k,i) - d_model_array(k,j)) / d_model_array(k,j))^2 / N);
                    end
                    
                    if (isnan(alpha_model_array(k,i)) || isnan(alpha_model_array(k,j)))
                    else
                        alpha_diff(i,j) = sqrt(alpha_diff(i,j)^2 + ((alpha_model_array(k,i) - alpha_model_array(k,j)) / alpha_model_array(k,j))^2 / N);
                    end
                end
            end
        end
    end
    
end


%% Weight gain least square criteria
%%{
if strcmp(alloy, 'Zry4') && T0 > T_transf
    %{
    % Comparison of experimental data
    
    % UW exp vs CP and C4 models
    idx = zeros(1, length(wg_exp_UW));
    wg_diff_UW_CP = 0;
    wg_diff_UW_C4 = 0;
    
    for k = 1:length(wg_exp_UW)
        t_curr = time_exp_UW(k);
        idx(k) = find(time == t_curr);
        wg_diff_UW_CP = wg_diff_UW_CP + ((wg_exp_UW(k) - wg_CP(idx(k))) / wg_CP(idx(k)))^2 / length(wg_exp_UW);
        wg_diff_UW_C4 = wg_diff_UW_C4 + ((wg_exp_UW(k) - wg(idx(k))) / wg_CP(idx(k)))^2 / length(wg_exp_UW);
    end
    %}
    
    %{
    % CP exp vs CP and C4 models
    idx = zeros(1, length(wg_exp_CP));
    wg_diff_CP_CP = 0;
    wg_diff_CP_C4 = 0;
    
    for k = 1:length(wg_exp_CP)
        t_curr = time_exp_CP(k);
        idx(k) = find(time == t_curr);
        wg_diff_CP_CP = wg_diff_CP_CP + ((wg_exp_CP(k) - wg_CP(idx(k))) / wg_CP(idx(k)))^2 / length(wg_exp_CP);
        wg_diff_CP_C4 = wg_diff_CP_C4 + ((wg_exp_CP(k) - wg(idx(k))) / wg_CP(idx(k)))^2 / length(wg_exp_CP);
    end
    %}
    
    % C4 model vs CP exp
    idx = zeros(1, length(wg_exp_CP));
    wg_diff_C4_CP_exp = 0;
    
    for k = 1:length(wg_exp_CP)
        t_curr = time_exp_CP(k);
        idx(k) = find(time == t_curr);
        wg_diff_C4_CP_exp = sqrt(wg_diff_C4_CP_exp^2 + ((wg(idx(k)) - wg_exp_CP(k)) / wg_exp_CP(k))^2 / length(wg_exp_CP));
    end
    
    % C4 model vs CP + UW exp
    if exist('wg_exp_UW', 'var') == 1
        idx = zeros(1, length(wg_exp_CP)+length(wg_exp_UW));
        wg_diff_C4_CPUW_exp = 0;
        
        for k = 1:length(wg_exp_CP)
            t_curr = time_exp_CP(k);
            idx(k) = find(time == t_curr);
            wg_diff_C4_CPUW_exp = sqrt(wg_diff_C4_CPUW_exp^2 + ((wg(idx(k)) - wg_exp_CP(k)) / wg_exp_CP(k))^2 / (length(wg_exp_CP)+length(wg_exp_UW)));
        end
        for k = 1:length(wg_exp_UW)
            t_curr = time_exp_UW(k);
            idx(k) = find(time == t_curr);
            wg_diff_C4_CPUW_exp = sqrt(wg_diff_C4_CPUW_exp^2 + ((wg(idx(k)) - wg_exp_UW(k)) / wg_exp_UW(k))^2 / (length(wg_exp_CP)+length(wg_exp_UW)));
        end
    end
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('Weight gain experimental comparison \n');
%     fprintf('Residual UW exp - CP model: %f\n', wg_diff_UW_CP);
%     fprintf('Residual UW exp - C4 model: %f\n', wg_diff_UW_C4);
%     fprintf('Residual CP exp - CP model: %f\n', wg_diff_CP_CP);
%     fprintf('Residual CP exp - C4 model: %f\n', wg_diff_CP_C4);
    fprintf('Residual C4 model - CP exp: %f\n', wg_diff_C4_CP_exp);
    if exist('wg_exp_UW', 'var') == 1
        fprintf('Residual C4 model - CP-UW exp: %f\n', wg_diff_C4_CPUW_exp);
    end
    %}
    
    % Comparison of models
    
    idx = zeros(1, length(wg_exp_CP));
    wg_diff_C4_CP_exp = 0;
    
    for k = 1:length(wg_exp_CP)
        t_curr = time_exp_CP(k);
        idx(k) = find(time == t_curr);
        wg_diff_C4_CP_exp = sqrt(wg_diff_C4_CP_exp^2 + ((wg(idx(k)) - wg_exp_CP(k)) / wg_exp_CP(k))^2 / length(wg_exp_CP));
    end
    
    wg_diff_BJ_CP = 0;
    wg_diff_Leist_CP = 0;
    wg_diff_WPI_CP = 0;
    wg_diff_Kawasaki_CP = 0;
    wg_diff_C4_CP = 0;
    
    wg_diff_BJ_C4 = 0;
    wg_diff_Leist_C4 = 0;
    wg_diff_WPI_C4 = 0;
    wg_diff_Kawasaki_C4 = 0;
    wg_diff_CP_C4 = 0;
    
    for k = 2:N+1
%         wg_diff_BJ_CP = sqrt(wg_diff_BJ_CP^2 + ((wg_BJ(k) - wg_CP_extended(k)) / wg_CP_extended(k))^2 / N);
%         wg_diff_Leist_CP = sqrt(wg_diff_Leist_CP^2 + ((wg_Leistikow(k) - wg_CP_extended(k)) / wg_CP_extended(k))^2 / N);
%         wg_diff_WPI_CP = sqrt(wg_diff_WPI_CP^2 + ((wg_WPI(k) - wg_CP_extended(k)) / wg_CP_extended(k))^2 / N);
%         wg_diff_Kawasaki_CP = sqrt(wg_diff_Kawasaki_CP^2 + ((wg_Kawasaki(k) - wg_CP_extended(k)) / wg_CP_extended(k))^2 / N);
        wg_diff_C4_CP = sqrt(wg_diff_C4_CP^2 + ((wg(k) - wg_CP_extended(k)) / wg_CP_extended(k))^2 / N);
        
%         wg_diff_BJ_C4 = sqrt(wg_diff_BJ_C4^2 + ((wg_BJ(k) - wg(k)) / wg(k))^2 / N);
%         wg_diff_Leist_C4 = sqrt(wg_diff_Leist_C4^2 + ((wg_Leistikow(k) - wg(k)) / wg(k))^2 / N);
%         wg_diff_WPI_C4 = sqrt(wg_diff_WPI_C4^2 + ((wg_WPI(k) - wg(k)) / wg(k))^2 / N);
%         wg_diff_Kawasaki_C4 = sqrt(wg_diff_Kawasaki_C4^2 + ((wg_Kawasaki(k) - wg(k)) / wg(k))^2 / N);
        wg_diff_CP_C4 = sqrt(wg_diff_CP_C4^2 + ((wg_CP_extended(k) - wg(k)) / wg(k))^2 / N);
    end
    
    fprintf('\nWeight gain model comparison \n');
%     fprintf('Residual Baker-Just model - CP model: %f\n', wg_diff_BJ_CP);
%     fprintf('Residual Leistikow model - CP model: %f\n', wg_diff_Leist_CP);
%     fprintf('Residual Biederman (WPI) model - CP model: %f\n', wg_diff_WPI_CP);
%     fprintf('Residual Kawasaki model - CP model: %f\n', wg_diff_Kawasaki_CP);
    fprintf('Residual C4 model - CP model: %f\n', wg_diff_C4_CP);
    
%     fprintf('\nResidual Baker-Just model - C4 model: %f\n', wg_diff_BJ_C4);
%     fprintf('Residual Leistikow model - C4 model: %f\n', wg_diff_Leist_C4);
%     fprintf('Residual Biederman (WPI) model - C4 model: %f\n', wg_diff_WPI_C4);
%     fprintf('Residual Kawasaki model - C4 model: %f\n', wg_diff_Kawasaki_C4);
%     fprintf('Residual CP model - C4 model: %f\n', wg_diff_CP_C4);
    
end
%}


%% Oxide thickness least square criteria
%%{

if strcmp(temperature, '360C')
    % Comparison of experimental data
    
    % C4 model vs exp
    idx = zeros(1, length(d_exp));
    d_diff_C4_exp = 0;
    
    for k = 1:length(d_exp)
        t_curr = 86400 * day_wg_exp(k);
        idx(k) = find(time == t_curr);
        d_diff_C4_exp = sqrt(d_diff_C4_exp^2 + ((1e4*d_total(idx(k)) - d_exp(k)) / d_exp(k))^2 / length(d_exp));
    end
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('Oxide thickness experimental comparison \n');
    fprintf('Residual C4 model - exp: %f\n', d_diff_C4_exp);
end

if strcmp(alloy, 'Zry4')
    %{
    % Comparison of experimental data
    
    % CP exp vs CP and C4 model
    idx = zeros(1, length(d_exp_CP));
    d_diff_CP_CP = 0;
    d_diff_CP_C4 = 0;
    
    for k = 1:length(d_exp_CP)
        t_curr = time_exp_CP(k);
        idx(k) = find(time == t_curr);
        d_diff_CP_CP = d_diff_CP_CP + ((d_exp_CP(k) - 1e4*d_CP(idx(k))) / (1e4*d_CP(idx(k))))^2 / length(d_exp_CP);
        d_diff_CP_C4 = d_diff_CP_C4 + ((d_exp_CP(k) - 1e4*d_protect(idx(k))) / (1e4*d_CP(idx(k))))^2 / length(d_exp_CP);
    end
    %}
    
    % C4 model vs CP exp
    idx = zeros(1, length(d_exp_CP));
    d_diff_C4_CP_exp = 0;
    
    for k = 1:length(d_exp_CP)
        t_curr = time_exp_CP(k);
        idx(k) = find(time == t_curr);
        d_diff_C4_CP_exp = sqrt(d_diff_C4_CP_exp^2 + ((1e4*d_protect(idx(k)) - d_exp_CP(k)) / d_exp_CP(k))^2 / length(d_exp_CP));
    end
    
    % C4 model vs CP + UW exp
    if exist('d_exp_UW', 'var') == 1
        idx = zeros(1, length(d_exp_CP)+length(d_exp_UW));
        d_diff_C4_CPUW_exp = 0;
        
        for k = 1:length(d_exp_CP)
            t_curr = time_exp_CP(k);
            idx(k) = find(time == t_curr);
            d_diff_C4_CPUW_exp = sqrt(d_diff_C4_CPUW_exp^2 + ((1e4*d_protect(idx(k)) - d_exp_CP(k)) / d_exp_CP(k))^2 / (length(d_exp_CP)+length(d_exp_UW)));
        end
        for k = 1:length(d_exp_UW)
            t_curr = time_exp_UW(k);
            idx(k) = find(time == t_curr);
            d_diff_C4_CPUW_exp = sqrt(d_diff_C4_CPUW_exp^2 + ((1e4*d_protect(idx(k)) - d_exp_UW(k)) / d_exp_UW(k))^2 / (length(d_exp_CP)+length(d_exp_UW)));
        end
    end
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('Oxide thickness experimental comparison \n');
%     fprintf('Residual CP exp - CP model: %f\n', d_diff_CP_CP);
%     fprintf('Residual CP exp - C4 model: %f\n', d_diff_CP_C4);
    fprintf('Residual C4 model - CP exp: %f\n', d_diff_C4_CP_exp);
    if exist('wg_exp_UW', 'var') == 1
        fprintf('Residual C4 model - CP-UW exp: %f\n', d_diff_C4_CPUW_exp);
    end
    %}
    
    % Comparison of models
    
    d_diff_BJ_CP = 0;
    d_diff_Leist_CP = 0;
    d_diff_WPI_CP = 0;
    d_diff_Kawasaki_CP = 0;
    d_diff_C4_CP = 0;
    
    d_diff_BJ_C4 = 0;
    d_diff_Leist_C4 = 0;
    d_diff_WPI_C4 = 0;
    d_diff_Kawasaki_C4 = 0;
    d_diff_C4_exp = 0;
    
    for k = 2:N+1
        d_diff_BJ_CP = sqrt(d_diff_BJ_CP^2 + ((1e4*d_BJ(k) - 1e4*d_CP_extended(k)) / (1e4*d_CP_extended(k)))^2 / N);
        d_diff_Leist_CP = sqrt(d_diff_Leist_CP^2 + ((1e4*d_Leistikow(k) - 1e4*d_CP_extended(k)) / (1e4*d_CP_extended(k)))^2 / N);
%         d_diff_WPI_CP = sqrt(d_diff_WPI_CP^2 + ((d_WPI(k) - 1e4*d_CP_extended(k)) / (1e4*d_CP_extended(k)))^2 / N);
        d_diff_Kawasaki_CP = sqrt(d_diff_Kawasaki_CP^2 + ((1e4*d_Kawasaki(k) - 1e4*d_CP_extended(k)) / (1e4*d_CP_extended(k)))^2 / N);
        d_diff_C4_CP = sqrt(d_diff_C4_CP^2 + ((1e4*d_protect(k) - 1e4*d_CP_extended(k)) / (1e4*d_CP_extended(k)))^2 / N);
        
        d_diff_BJ_C4 = sqrt(d_diff_BJ_C4^2 + ((1e4*d_BJ(k) - 1e4*d_protect(k)) / (1e4*d_protect(k)))^2 / N);
        d_diff_Leist_C4 = sqrt(d_diff_Leist_C4^2 + ((1e4*d_Leistikow(k) - 1e4*d_protect(k)) / (1e4*d_protect(k)))^2 / N);
%         d_diff_WPI_C4 = sqrt(d_diff_WPI_C4^2 + ((d_WPI(k) - 1e4*d_protect(k)) / (1e4*d_protect(k)))^2 / N);
        d_diff_Kawasaki_C4 = sqrt(d_diff_Kawasaki_C4^2 + ((1e4*d_Kawasaki(k) - 1e4*d_protect(k)) / (1e4*d_protect(k)))^2 / N);
        d_diff_C4_exp = sqrt(d_diff_C4_exp^2 + ((1e4*d_CP_extended(k) - 1e4*d_protect(k)) / (1e4*d_protect(k)))^2 / N);
    end
    
    fprintf('\nOxide thickness model comparison \n');
%     fprintf('Residual Baker-Just model - CP model: %f\n', d_diff_BJ_CP);
%     fprintf('Residual Leistikow model - CP model: %f\n', d_diff_Leist_CP);
%     fprintf('Residual Biederman (WPI) model - CP model: %f\n', d_diff_WPI_CP);
%     fprintf('Residual Kawasaki model - CP model: %f\n', d_diff_Kawasaki_CP);
    fprintf('Residual C4 model - CP model: %f\n', d_diff_C4_CP);
    
%     fprintf('\nResidual Baker-Just model - C4 model: %f\n', d_diff_BJ_C4);
%     fprintf('Residual Leistikow model - C4 model: %f\n', d_diff_Leist_C4);
%     fprintf('Residual Biederman (WPI) model - C4 model: %f\n', d_diff_WPI_C4);
%     fprintf('Residual Kawasaki model - C4 model: %f\n', d_diff_Kawasaki_C4);
    fprintf('Residual CP model - C4 model: %f\n', d_diff_C4_exp);
end
%}


%% Alpha thickness least square criteria
%%{
if strcmp(alloy, 'Zry4')
    %{
    % Comparison of experimental data
    
    % C4 model vs CP exp
    idx = zeros(1, length(wg_exp_CP));
    wg_diff_C4_CPexp = 0;
    
    for k = 1:length(wg_exp_CP)
        t_curr = time_exp_CP(k);
        idx(k) = find(time == t_curr);
        wg_diff_C4_CPexp = sqrt(wg_diff_CP_C4^2 + ((wg(idx(k)) - wg_exp_CP(k)) / wg_exp_CP(k))^2 / length(wg_exp_CP));
    end
    %}
    
    % C4 model vs CP exp
    idx = zeros(1, length(alpha_exp_CP));
    alpha_diff_C4_CP_exp = 0;
    
    for k = 1:length(alpha_exp_CP)
        t_curr = time_exp_CP(k);
        idx(k) = find(time == t_curr);
        alpha_diff_C4_CP_exp = sqrt(alpha_diff_C4_CP_exp^2 + ((1e4*alpha(idx(k)) - alpha_exp_CP(k)) / alpha_exp_CP(k))^2 / length(alpha_exp_CP));
    end
    
    % C4 model vs CP + UW exp
    if exist('alpha_exp_UW', 'var') == 1
        idx = zeros(1, length(alpha_exp_CP)+length(alpha_exp_UW));
        alpha_diff_C4_CPUW_exp = 0;
        
        for k = 1:length(alpha_exp_CP)
            t_curr = time_exp_CP(k);
            idx(k) = find(time == t_curr);
            alpha_diff_C4_CPUW_exp = sqrt(alpha_diff_C4_CPUW_exp^2 + ((1e4*alpha(idx(k)) - alpha_exp_CP(k)) / alpha_exp_CP(k))^2 / (length(alpha_exp_CP)+length(alpha_exp_UW)));
        end
        for k = 1:length(alpha_exp_UW)
            t_curr = time_exp_UW(k);
            idx(k) = find(time == t_curr);
            alpha_diff_C4_CPUW_exp = sqrt(alpha_diff_C4_CPUW_exp^2 + ((1e4*alpha(idx(k)) - alpha_exp_UW(k)) / alpha_exp_UW(k))^2 / (length(alpha_exp_CP)+length(alpha_exp_UW)));
        end
    end
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('Alpha thickness experimental comparison \n');
%     fprintf('Residual CP exp - CP model: %f\n', alpha_diff_CP_CP);
%     fprintf('Residual CP exp - C4 model: %f\n', alpha_diff_CP_C4);
    fprintf('Residual C4 model - CP exp: %f\n', alpha_diff_C4_CP_exp);
    if exist('wg_exp_UW', 'var') == 1
        fprintf('Residual C4 model - CP-UW exp: %f\n', alpha_diff_C4_CPUW_exp);
    end
    %}
    
    % Comparison of models
    
    alpha_diff_WPI_CP = 0;
    alpha_diff_Kawasaki_CP = 0;
    alpha_diff_C4_CP = 0;
    
    alpha_diff_WPI_C4 = 0;
    alpha_diff_Kawasaki_C4 = 0;
    alpha_diff_CP_C4 = 0;
    
    for k = 3:N+1
%         alpha_diff_WPI_CP = sqrt(alpha_diff_WPI_CP^2 + ((alpha_WPI(k) - 1e4*alpha_CP_extended(k)) / (1e4*alpha_CP_extended(k)))^2 / (N-1));
        alpha_diff_Kawasaki_CP = sqrt(alpha_diff_Kawasaki_CP^2 + ((1e4*alpha_Kawasaki(k) - 1e4*alpha_CP_extended(k)) / (1e4*alpha_CP_extended(k)))^2 / (N-1));
        alpha_diff_C4_CP = sqrt(alpha_diff_C4_CP^2 + ((1e4*alpha(k) - 1e4*alpha_CP_extended(k)) / (1e4*alpha_CP_extended(k)))^2 / (N-1));
        
%         alpha_diff_WPI_C4 = sqrt(alpha_diff_WPI_C4^2 + ((alpha_WPI(k) - 1e4*alpha(k)) / (1e4*alpha(k)))^2 / (N-1));
        alpha_diff_Kawasaki_C4 = sqrt(alpha_diff_Kawasaki_C4^2 + ((1e4*alpha_Kawasaki(k) - 1e4*alpha(k)) / (1e4*alpha(k)))^2 / (N-1));
        alpha_diff_CP_C4 = sqrt(alpha_diff_CP_C4^2 + ((1e4*alpha_CP_extended(k) - 1e4*alpha(k)) / (1e4*alpha(k)))^2 / (N-1));
    end
    
    fprintf('\nAlpha thickness model comparison \n');
%     fprintf('Residual Biederman (WPI) model - CP model: %f\n', alpha_diff_WPI_CP);
%     fprintf('Residual Kawasaki model - CP model: %f\n', alpha_diff_Kawasaki_CP);
    fprintf('Residual C4 model - CP model: %f\n', alpha_diff_C4_CP);
    
%     fprintf('\nResidual Biederman (WPI) model - C4 model: %f\n', alpha_diff_WPI_C4);
%     fprintf('Residual Kawasaki model - C4 model: %f\n', alpha_diff_Kawasaki_C4);
    fprintf('Residual CP model - C4 model: %f\n', alpha_diff_CP_C4);
end
%}