%% C4_MOOSE_comparison
% Compute gradients at the interface, interfaces position, NRMSE on oxide thcikness and alpha thickness
% and alpha/beta interface velocity for Matlab and MOOSE isotherm HT oxidation data. 

%%Author : Louis Bailly-Salins
%%Email : baillysalins@wisc.edu

%%Last updated : 08/05/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gradients at the interfaces (mixed with old plots for 1200C)

%Retrieving data from MOOSE postprocessor csv file (IC5 bis (cst between 2 alpha additional interfaces), dx=3µm, new a/b velocity)
% MOOSE_data_new = csvread('MOOSE_data/2interfaces_Da10_Db60_IC5bis_3um_newvelocity.csv',1,0);
% time_MOOSE_500 = MOOSE_data(1:480,1);
% grad_ab_MOOSE_new = 10000*MOOSE_data_new(1:480,2);
% grad_aox_MOOSE_new = 10000*MOOSE_data_new(1:480,3);
% grad_ba_MOOSE_new = 10000*MOOSE_data_new(1:480,4);
% % Transforms data a bit : shift 1s earlier (because executed on timestep begin, before solving)
% time_MOOSE_500 = time_MOOSE_500 - 1;
% % Erase first value (not right)
% time_MOOSE_500(1) = NaN;
% grad_ab_MOOSE_new(1) = NaN;
% grad_aox_MOOSE_new(1) = NaN;
% grad_ba_MOOSE_new(1) = NaN;

tmax = 1500;

if temperature == '1000C'
    MOOSE_data = csvread('MOOSE_data/T_check/T_1273K.csv',1,0);
end

if temperature == '1100C'
    MOOSE_data = csvread('MOOSE_data/T_check/T_1373K.csv',1,0);
end
    

if temperature == '1200C' 
%     % Retrieving data from MOOSE postprocessor csv file (IC5 bis (cst between 2 alpha additional interfaces), dx=1µm)
%     MOOSE_data_IC5b = csvread('MOOSE_data/right_vab/2interfaces_Da10_Db60_IC5bis_1um_newvelocity.csv',1,0);
%     time_MOOSE_500 = MOOSE_data_IC5b(1:480,1);
%     grad_ab_MOOSE_IC5b = 10000*MOOSE_data_IC5b(1:480,2);
%     grad_aox_MOOSE_IC5b = 10000*MOOSE_data_IC5b(1:480,3);
%     grad_ba_MOOSE_IC5b = 10000*MOOSE_data_IC5b(1:480,4);
%     % % Transforms data a bit : shift 1s earlier (because executed on timestep begin, before solving)
%     time_MOOSE_500 = time_MOOSE_500 - 1;
%     % Erase first value (not right)
%     time_MOOSE_500(1) = NaN;
%     grad_ab_MOOSE_IC5b(1) = NaN;
%     grad_aox_MOOSE_IC5b(1) = NaN;
%     grad_ba_MOOSE_IC5b(1) = NaN;
% 
%     % Retrieving data from MOOSE postprocessor csv file (IC4, dx=1µm)
%     MOOSE_data_IC4 = csvread('MOOSE_data/right_vab/2interfaces_Da10_Db60_IC4_1um_newvelocity.csv',1,0);
%     grad_ab_MOOSE_IC4 = 10000*MOOSE_data_IC4(1:480,2);
%     grad_aox_MOOSE_IC4 = 10000*MOOSE_data_IC4(1:480,3);
%     grad_ba_MOOSE_IC4 = 10000*MOOSE_data_IC4(1:480,4);
%     % Erase first value (not right)
%     grad_ab_MOOSE_IC4(1) = NaN;
%     grad_aox_MOOSE_IC4(1) = NaN;
%     grad_ba_MOOSE_IC4(1) = NaN;
% 
%     % Retrieving data from MOOSE postprocessor csv file (ICconst, dx=1µm)
%     MOOSE_data_ICc = csvread('MOOSE_data/right_vab/moving_oxide_C4_weak_ICconst_out.csv',1,0);
%     grad_ab_MOOSE_ICc = 10000*MOOSE_data_ICc(1:480,2);
%     grad_aox_MOOSE_ICc = 10000*MOOSE_data_ICc(1:480,3);
%     grad_ba_MOOSE_ICc = 10000*MOOSE_data_ICc(1:480,4);
%     % Erase first value (not right)
%     grad_ab_MOOSE_ICc(1) = NaN;
%     grad_aox_MOOSE_ICc(1) = NaN;
%     grad_ba_MOOSE_ICc(1) = NaN;

%     % Retrieving data from MOOSE postprocessor csv file (IC5, dx=1µm)
%     MOOSE_data_IC5 = csvread('MOOSE_data/right_vab/moving_oxide_C4_weak_IC5_out.csv',1,0);
%     grad_ab_MOOSE_IC5 = 10000*MOOSE_data_IC5(1:480,2);
%     grad_aox_MOOSE_IC5 = 10000*MOOSE_data_IC5(1:480,3);
%     grad_ba_MOOSE_IC5 = 10000*MOOSE_data_IC5(1:480,4);
%     %   Erase first value (not right)
%     grad_ab_MOOSE_IC5(1) = NaN;
%     grad_aox_MOOSE_IC5(1) = NaN;
%     grad_ba_MOOSE_IC5(1) = NaN;
% 
%     % Retrieving data from MOOSE postprocessor csv file (ICconstlin, dx=1µm)
%     MOOSE_data_ICcl = csvread('MOOSE_data/right_vab/moving_oxide_C4_weak_ICconstlin_out.csv',1,0);
%     grad_ab_MOOSE_ICcl = 10000*MOOSE_data_ICcl(1:480,2);
%     grad_aox_MOOSE_ICcl = 10000*MOOSE_data_ICcl(1:480,3);
%     grad_ba_MOOSE_ICcl = 10000*MOOSE_data_ICcl(1:480,4);
%     % Erase first value (not right)
%     grad_ab_MOOSE_ICcl(1) = NaN;
%     grad_aox_MOOSE_ICcl(1) = NaN;
%     grad_ba_MOOSE_ICcl(1) = NaN;
% 
%     % Retrieving data from MOOSE postprocessor csv file (IClin, dx=1µm)
%     MOOSE_data_ICl = csvread('MOOSE_data/right_vab/moving_oxide_C4_weak_IClin_out.csv',1,0); %ICl for IClin (l the letter, not 1 the number)
%     grad_ab_MOOSE_ICl = 10000*MOOSE_data_ICl(1:480,2);
%     grad_aox_MOOSE_ICl = 10000*MOOSE_data_ICl(1:480,3);
%     grad_ba_MOOSE_ICl = 10000*MOOSE_data_ICl(1:480,4);
%     % Erase first value (not right)
%     grad_ab_MOOSE_ICl(1) = NaN;
%     grad_aox_MOOSE_ICl(1) = NaN;
%     grad_ba_MOOSE_ICl(1) = NaN;
% 
%     % Retrieving data from MOOSE postprocessor csv file (ICconst, start time = 10s, dx=1µm)
%     MOOSE_data_t10 = csvread('MOOSE_data/right_vab/moving_oxide_C4_weak_ICconst_10s_out.csv',1,0);
%     time_t10 = MOOSE_data_t10(1:490,1);
%     time_t10 = time_t10 - 1;
%     grad_ab_MOOSE_t10 = 10000*MOOSE_data_t10(1:490,2);
%     grad_aox_MOOSE_t10 = 10000*MOOSE_data_t10(1:490,3);
%     grad_ba_MOOSE_t10 = 10000*MOOSE_data_t10(1:490,4);
%     % Erase first value (not right)
%     time_t10(1) = NaN;
%     grad_ab_MOOSE_t10(1) = NaN;
%     grad_aox_MOOSE_t10(1) = NaN;
%     grad_ba_MOOSE_t10(1) = NaN;

%     % Retrieving data from MOOSE postprocessor csv file (ICconst, start time = 30s, dx=1µm)
%     MOOSE_data_t30 = csvread('MOOSE_data/right_vab/moving_oxide_C4_weak_ICconst_30s_out.csv',1,0);
%     time_t30 = MOOSE_data_t30(1:470,1);
%     time_t30 = time_t30 - 1;
%     grad_ab_MOOSE_t30 = 10000*MOOSE_data_t30(1:470,2);
%     grad_aox_MOOSE_t30 = 10000*MOOSE_data_t30(1:470,3);
%     grad_ba_MOOSE_t30 = 10000*MOOSE_data_t30(1:470,4);
%     % Erase first value (not right)
%     time_t30(1) = NaN;
%     grad_ab_MOOSE_t30(1) = NaN;
%     grad_aox_MOOSE_t30(1) = NaN;
%     grad_ba_MOOSE_t30(1) = NaN;

    % plot(time_MOOSE_500, grad_aox_MOOSE_new,'LineStyle', '--', 'Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=3\mum, IC5bis, new velocity)');
    % plot(time_MOOSE_500, grad_ab_MOOSE_new,'LineStyle', '--', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=3\mum, IC5bis, new velocity)');
    % plot(time_MOOSE_500, grad_ba_MOOSE_new,'LineStyle', '--', 'Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=3\mum, IC5bis, new velocity)');

    % plot(time_MOOSE_500, grad_aox_MOOSE_IC5b,'Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IC5bis)');
    % plot(time_MOOSE_500, grad_ab_MOOSE_IC5b, 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IC5bis)');
    % plot(time_MOOSE_500, grad_ba_MOOSE_IC5b,'Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IC5bis)');

    % plot(time_MOOSE_500, grad_aox_MOOSE_IC4,'LineStyle', '--','Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IC4)');
    % plot(time_MOOSE_500, grad_ab_MOOSE_IC4,'LineStyle', '--', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IC4)');
    % plot(time_MOOSE_500, grad_ba_MOOSE_IC4,'LineStyle', '--','Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IC4)');

    % plot(time_MOOSE_500, grad_aox_MOOSE_IC5,'LineStyle', ':','Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IC5)');
    % plot(time_MOOSE_500, grad_ab_MOOSE_IC5,'LineStyle', ':', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IC5)');
    % plot(time_MOOSE_500, grad_ba_MOOSE_IC5,'LineStyle', ':','Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IC5)');

    % plot(time_MOOSE_500, grad_aox_MOOSE_ICl,'LineStyle', '-.','Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IClin)');
    % plot(time_MOOSE_500, grad_ab_MOOSE_ICl,'LineStyle', '-.', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IClin)');
    % plot(time_MOOSE_500, grad_ba_MOOSE_ICl,'LineStyle', '-.','Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IClin)');

    % plot(time_MOOSE_500, grad_aox_MOOSE_ICcl,'LineStyle', 'None','Marker','x','Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, ICconstlin)');
    % plot(time_MOOSE_500, grad_ab_MOOSE_ICcl,'LineStyle', 'None','Marker','x', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, ICconstlin)');
    % plot(time_MOOSE_500, grad_ba_MOOSE_ICcl,'LineStyle', 'None','Marker','x','Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, ICconstlin)');

%     plot(time_t10,grad_aox_MOOSE_t10,'LineStyle', '--','Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, ICconst, t_{init}=10s)');
%     plot(time_t10, grad_ab_MOOSE_t10,'LineStyle', '--', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, ICconst, t_{init}=10s)');
%     plot(time_t10,grad_ba_MOOSE_t10 ,'LineStyle', '--','Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, ICconst, t_{init}=10s)');
% 
 
%     plot(time_t30,grad_aox_MOOSE_t30,'LineStyle', ':','Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, ICconst, t_{init}=30s)');
%     plot(time_t30, grad_ab_MOOSE_t30,'LineStyle', ':', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, ICconst, t_{init}=30s)');
%     plot(time_t30,grad_ba_MOOSE_t30 ,'LineStyle', ':','Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, ICconst, t_{init}=30s)');

    % Retrieving data from MOOSE postprocessor csv file (ICconst, dx=1µm)
    MOOSE_data = csvread('MOOSE_data/T_check/T_1473K.csv',1,0);
end

if temperature == '1300C' 
    tmax = 200;
    % Retrieving data from MOOSE postprocessor csv file
    MOOSE_data = csvread('MOOSE_data/T_check/T_1573K.csv',1,0);
end

if temperature == '1400C' 
    tmax = 100
    % Retrieving data from MOOSE postprocessor csv file
    MOOSE_data = csvread('MOOSE_data/T_check/T_1673K.csv',1,0);   
end

if temperature == '1500C' 
    tmax = 60;
    % Retrieving data from MOOSE postprocessor csv file
    MOOSE_data = csvread('MOOSE_data/T_check/T_1773K.csv',1,0);
end

time_MOOSE = MOOSE_data(1:(tmax-20+1),1);
grad_ab_MOOSE = 10000*MOOSE_data(1:(tmax-20+1),2);
grad_aox_MOOSE = 10000*MOOSE_data(1:(tmax-20+1),3);
grad_ba_MOOSE = 10000*MOOSE_data(1:(tmax-20+1),4);
grad_ab_MOOSE(1) = NaN;
grad_aox_MOOSE(1) = NaN;
grad_ba_MOOSE(1) = NaN;

%% MOOSE Interfaces position 

xi_ox_MOOSE = MOOSE_data(1:(tmax-20+1),6);
xi_ab_MOOSE = MOOSE_data(1:(tmax-20+1),5);

%% NRMSE on oxide and alpha layer thickness, MOOSE vs Matlab

d_MOOSE = 0.0001 * PBR * abs(600-xi_ox_MOOSE);     % MOOSE oxide thickness [cm]
alpha_MOOSE = 0.0001 * abs(xi_ox_MOOSE-xi_ab_MOOSE);  % MOOSE alpha thickness [cm]

N_MOOSE = length(xi_ox_MOOSE);
d_diff_MOOSE = zeros(1,N_MOOSE);
alpha_diff_MOOSE = zeros(1,N_MOOSE);
instant_d_diff = zeros(1,N_MOOSE);
instant_alpha_diff = zeros(1,N_MOOSE);

instant_d_diff(1) = abs((d_MOOSE(1)-d_total(201))/d_total(201));
instant_alpha_diff(1) = abs((alpha_MOOSE(1)-alpha(201))/alpha(201));
d_diff_MOOSE(1) = instant_d_diff(1) ;
alpha_diff_MOOSE(1) = instant_alpha_diff(1);

for k=2:N_MOOSE;
   k_matlab = (k+19)*10+1;
   instant_d_diff(k) = abs((d_MOOSE(k)-d_total(k_matlab))/d_total(k_matlab));
   instant_alpha_diff(k) = abs((alpha_MOOSE(k)-alpha(k_matlab))/alpha(k_matlab));
   %instant_d_diff(k) = abs((d_MOOSE(k)-d_total(k))/d_total(k));   %for
   %dt=1 in Matlab
   %instant_alpha_diff(k) = abs((alpha_MOOSE(k)-alpha(k))/alpha(k)); %for
   %dt=1 in Matlab
   d_diff_MOOSE(k) = sqrt((k-1)/k*d_diff_MOOSE(k-1)^2 + instant_d_diff(k)^2/k);
   alpha_diff_MOOSE(k) = sqrt((k-1)/k*alpha_diff_MOOSE(k-1)^2 + instant_alpha_diff(k)^2/k);
end


%% Alpha/beta interface velocity

v_ab = (Db*grad_ba - Da*grad_ab)/(Co_a_b-Co_b_a);                       % Matlab alpha/beta velocity [cm/s]
v_ab_MOOSE = Czr*(Db*grad_ba_MOOSE - Da*grad_ab_MOOSE)/(Co_a_b-Co_b_a); % MOOSE alpha/beta velocity [cm/s]
v_ab_MOOSE(1)=NaN;

v_ab_obs = zeros(1,N);
v_ab_obs(1) = NaN; %t=0
v_ab_obs(2) = NaN; %t=0.1s
v_ab_obs(3) = NaN; %t=0.2s
for k=5:N+1
    v_ab_obs(k-1)=(xi_ab(k)-xi_ab(k-1))/dt;
end


%% Vacancy flux as a function of delta
%{
d_MOOSE = [1.248 8.53161 8.99318 23.3466 30.0667 41.1152 49.9815 61.9487 74.2039 77.2495];
Jv_MOOSE = 10^21*[325.594 47.6277 45.1832 17.4047 13.5146 9.88299 8.12984 6.55931 5.47601 5.26012];
%}

%% Gradients at the interfaces with wrong alpha/beta velocity (wrong and mixed with plot)
%{
fig_grads = figure('Name', 'Gradients at the interfaces');

% Retrieving data from MOOSE postprocessor csv file (dx=3µm)
MOOSE_data = csvread('MOOSE_data/2_interfaces_Da10_Db60_IC4int.csv',1,0);
time_MOOSE = MOOSE_data(:,1);
grad_ab_MOOSE = 10000*MOOSE_data(:,2);
grad_aox_MOOSE = 10000*MOOSE_data(:,3);
grad_ba_MOOSE = 10000*MOOSE_data(:,4);
% Transforms data a bit : shift 1s earlier (because executed on timestep begin, before solving)
time_MOOSE = time_MOOSE - 1;
% Erase first value (not right)
time_MOOSE(1) = NaN;
%time_MOOSE(2) = NaN;
grad_ab_MOOSE(1) = NaN;
grad_aox_MOOSE(1) = NaN;
grad_ba_MOOSE(1) = NaN;
%grad_ab_MOOSE(2) = NaN;
%grad_aox_MOOSE(2) = NaN;
%grad_ba_MOOSE(2) = NaN;

% Retrieving data from MOOSE postprocessor csv file (dx=1µm)
MOOSE_data1 = csvread('MOOSE_data/2_interfaces_Da10_Db60_IC4_1um.csv',1,0);
grad_ab_MOOSE1 = 10000*MOOSE_data1(:,2);
grad_aox_MOOSE1 = 10000*MOOSE_data1(:,3);
grad_ba_MOOSE1 = 10000*MOOSE_data1(:,4);
% Erase first value (not right)
grad_ab_MOOSE1(1) = NaN;
grad_aox_MOOSE1(1) = NaN;
grad_ba_MOOSE1(1) = NaN;
%grad_ab_MOOSE1(2) = NaN;
%grad_aox_MOOSE1(2) = NaN;
%grad_ba_MOOSE1(2) = NaN;

% Retrieving data from MOOSE postprocessor csv file (IC5, dx=1µm)
MOOSE_data_IC5 = csvread('MOOSE_data/2interfaces_Da10_Db60_IC5_1um.csv',1,0);
grad_ab_MOOSE_IC5 = 10000*MOOSE_data_IC5(:,2);
grad_aox_MOOSE_IC5 = 10000*MOOSE_data_IC5(:,3);
grad_ba_MOOSE_IC5 = 10000*MOOSE_data_IC5(:,4);
% Erase first value (not right)
grad_ab_MOOSE_IC5(1) = NaN;
grad_aox_MOOSE_IC5(1) = NaN;
grad_ba_MOOSE_IC5(1) = NaN;
%grad_ab_MOOSE_IC5(2) = NaN;
%grad_aox_MOOSE_IC5(2) = NaN;
%grad_ba_MOOSE_IC5(2) = NaN;

% Retrieving data from MOOSE postprocessor csv file (IC5 bis (cst between 2 alpha additional interfaces), dx=1µm)
MOOSE_data_IC5b = csvread('MOOSE_data/2interfaces_Da10_Db60_IC5bis_1um.csv',1,0);
time_MOOSE_200 = MOOSE_data(1:180,1);
grad_ab_MOOSE_IC5b = 10000*MOOSE_data_IC5b(1:180,2);
grad_aox_MOOSE_IC5b = 10000*MOOSE_data_IC5b(1:180,3);
grad_ba_MOOSE_IC5b = 10000*MOOSE_data_IC5b(1:180,4);
% Transforms data a bit : shift 1s earlier (because executed on timestep begin, before solving)
time_MOOSE_200 = time_MOOSE_200 - 1;
% Erase first value (not right)
time_MOOSE_200(1)= NaN;
grad_ab_MOOSE_IC5b(1) = NaN;
grad_aox_MOOSE_IC5b(1) = NaN;
grad_ba_MOOSE_IC5b(1) = NaN;

% Retrieving data from MOOSE postprocessor csv file (1 interface, 500s, dx=1um)
MOOSE_data_1interf = csvread('MOOSE_data/1_interface_Da10_Db60_1um_500s.csv',1,0);
time_1interf = MOOSE_data_1interf(:,1); % Have to retrieve the time as it's shorter
% Transforms data a bit : shift 1s earlier (because executed on timestep begin, before solving)
time_1interf = time_1interf - 1;
grad_aox_1interf = 10000*MOOSE_data_1interf(:,2);
% Erase first value (not right)
time_1interf(1) = NaN;
%time_1interf(2) = NaN;
grad_aox_1interf(1) = NaN;
%grad_aox_1interf(2) = NaN;

% Retrieving data from MOOSE postprocessor csv file (IC5 bis (cst between 2 alpha additional interfaces), dx=0.4µm)
MOOSE_data_04 = csvread('MOOSE_data/2interfaces_Da10_Db60_IC5bis_04um.csv',1,0);
time_MOOSE_100 = MOOSE_data(1:80,1);
grad_ab_MOOSE_04 = 10000*MOOSE_data_04(1:80,2);
grad_aox_MOOSE_04 = 10000*MOOSE_data_04(1:80,3);
grad_ba_MOOSE_04 = 10000*MOOSE_data_04(1:80,4);
% Transforms data a bit : shift 1s earlier (because executed on timestep begin, before solving)
time_MOOSE_100 = time_MOOSE_100 - 1;
% Erase first value (not right)
time_MOOSE_100(1)= NaN;
grad_ab_MOOSE_04(1) = NaN;
grad_aox_MOOSE_04(1) = NaN;
grad_ba_MOOSE_04(1) = NaN;

% Retrieving data from MOOSE postprocessor csv file (IC5 bis (cst between 2 alpha additional interfaces), dx=0.2µm)
MOOSE_data_02 = csvread('MOOSE_data/2interfaces_Da10_Db60_IC5bis_02um.csv',1,0);
grad_ab_MOOSE_02 = 10000*MOOSE_data_02(1:80,2);
grad_aox_MOOSE_02 = 10000*MOOSE_data_02(1:80,3);
grad_ba_MOOSE_02 = 10000*MOOSE_data_02(1:80,4);
% Erase first value (not right)
grad_ab_MOOSE_02(1) = NaN;
grad_aox_MOOSE_02(1) = NaN;
grad_ba_MOOSE_02(1) = NaN;

% Retrieving data from MOOSE postprocessor csv file (IC5 bis (cst between 2 alpha additional interfaces), dx=3µm, new a/b velocity)
MOOSE_data_new = csvread('MOOSE_data/2interfaces_Da10_Db60_IC5bis_3um_newvelocity.csv',1,0);
time_MOOSE_500 = MOOSE_data(1:480,1);
grad_ab_MOOSE_new = 10000*MOOSE_data_new(1:480,2);
grad_aox_MOOSE_new = 10000*MOOSE_data_new(1:480,3);
grad_ba_MOOSE_new = 10000*MOOSE_data_new(1:480,4);
% Transforms data a bit : shift 1s earlier (because executed on timestep begin, before solving)
time_MOOSE_500 = time_MOOSE_500 - 1;
% Erase first value (not right)
time_MOOSE_500(1) = NaN;
grad_ab_MOOSE_new(1) = NaN;
grad_aox_MOOSE_new(1) = NaN;
grad_ba_MOOSE_new(1) = NaN;

% Retrieving data from MOOSE postprocessor csv file (IC5 bis (cst between 2 alpha additional interfaces), dx=1µm, new a/b velocity)
MOOSE_data_new1 = csvread('MOOSE_data/2interfaces_Da10_Db60_IC5bis_1um_newvelocity.csv',1,0);
grad_ab_MOOSE_new1 = 10000*MOOSE_data_new1(1:480,2);
grad_aox_MOOSE_new1 = 10000*MOOSE_data_new1(1:480,3);
grad_ba_MOOSE_new1 = 10000*MOOSE_data_new1(1:480,4);
% Erase first value (not right)
grad_ab_MOOSE_new1(1) = NaN;
grad_aox_MOOSE_new1(1) = NaN;
grad_ba_MOOSE_new1(1) = NaN;


hold on;
plot(time_plot, grad_aox/Czr, 'Color', [1 0.5 0], 'DisplayName', 'Gradient at \alpha/oxide interface');
plot(time_plot, grad_ab/Czr, 'b', 'DisplayName', 'Gradient at \alpha/\beta interface');
plot(time_plot, grad_ba/Czr, 'Color', [0 0.6 0], 'DisplayName', 'Gradient at \beta/\alpha interface');

% plot(time_MOOSE, grad_aox_MOOSE, 'Color', [1 0.5 0], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface(\Deltax=3\mum)');
% plot(time_MOOSE, grad_ab_MOOSE,'b', 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface (\Deltax=3\mum)');
% plot(time_MOOSE, grad_ba_MOOSE, 'Color', [0 0.6 0], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface (\Deltax=3\mum)');

% plot(time_MOOSE, grad_aox_MOOSE1,'LineStyle', ':', 'Color', [1 0.5 0], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IC4)');
% plot(time_MOOSE, grad_ab_MOOSE1,'b:', 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IC4)');
% plot(time_MOOSE, grad_ba_MOOSE1,'LineStyle', ':', 'Color', [0 0.6 0], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IC4)');

%plot(time_1interf, grad_aox_1interf,'k-.', 'DisplayName', 'MOOSE Gradient at single \alpha/oxide interface(\Deltax=1\mum)');

% plot(time_MOOSE, grad_aox_MOOSE_IC5,'LineStyle', '--', 'Color', [1 0.5 0], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IC5)');
% plot(time_MOOSE, grad_ab_MOOSE_IC5,'b--', 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IC5)');
% plot(time_MOOSE, grad_ba_MOOSE_IC5,'LineStyle', '--', 'Color', [0 0.6 0], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IC5)');

plot(time_MOOSE_200, grad_aox_MOOSE_IC5b,'LineStyle', '--', 'Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IC5bis)');
plot(time_MOOSE_200, grad_ab_MOOSE_IC5b,'LineStyle', '--', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IC5bis)');
plot(time_MOOSE_200, grad_ba_MOOSE_IC5b,'LineStyle', '--', 'Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IC5bis)');

% plot(time_MOOSE_100, grad_aox_MOOSE_04,'LineStyle', ':', 'Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=0.4\mum, IC5bis)');
% plot(time_MOOSE_100, grad_ab_MOOSE_04,'LineStyle', ':', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=0.4\mum, IC5bis)');
% plot(time_MOOSE_100, grad_ba_MOOSE_04,'LineStyle', ':', 'Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=0.4\mum, IC5bis)');

% plot(time_MOOSE_100, grad_aox_MOOSE_02, 'Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=0.2\mum, IC5bis)');
% plot(time_MOOSE_100, grad_ab_MOOSE_02, 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=0.2\mum, IC5bis)');
% plot(time_MOOSE_100, grad_ba_MOOSE_02, 'Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=0.2\mum, IC5bis)');

% plot(time_MOOSE_500, grad_aox_MOOSE_new,'LineStyle', '--', 'Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=3\mum, IC5bis, new velocity)');
% plot(time_MOOSE_500, grad_ab_MOOSE_new,'LineStyle', '--', 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=3\mum, IC5bis, new velocity)');
% plot(time_MOOSE_500, grad_ba_MOOSE_new,'LineStyle', '--', 'Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=3\mum, IC5bis, new velocity)');

plot(time_MOOSE_500, grad_aox_MOOSE_new1,'Color', [0.929 0.694 0.125], 'DisplayName', 'MOOSE Gradient at \alpha/oxide interface (\Deltax=1\mum, IC5bis, new velocity)');
plot(time_MOOSE_500, grad_ab_MOOSE_new1, 'Color', [0.301 0.745 0.933], 'DisplayName', 'MOOSE Gradient at \alpha/\beta interface(\Deltax=1\mum, IC5bis, new velocity)');
plot(time_MOOSE_500, grad_ba_MOOSE_new1,'Color', [0.466 0.674 0.188], 'DisplayName', 'MOOSE Gradient at \beta/\alpha interface(\Deltax=1\mum, IC5bis, new velocity)');


xlim([0 500]);
xlabel(strcat('Exposure time [', time_unit, ']'));
ylim([0 300]);
ylabel('Gradients of Co/Czr at the interface [/cm]');
title('Gradients at the interfaces evolution with exposure time');
legend('-DynamicLegend', 'Location', 'NorthEast');
%}