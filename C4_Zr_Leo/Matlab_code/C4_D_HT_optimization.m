%% C4_D_HT_optimization
%%Find the best high temperature oxygen diffusion coefficients

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Dichotomy search

C4_parameter;

tolerance = 5e-2;
% Da_bound = [0.196*exp(-42500/(RR*T0)) 0.196*exp(-39500/(RR*T0))];
% Db_bound = [1.71e-2*exp(-29300./(RR*T0)) 0.0453*exp(-25800/(RR*T0))];
Da_bound = [1e-15 1];
Db_bound = [1e-15 1];
max_iter = 500;
Da = 0.196 * exp(-41000/(RR*T0));       % Mallett
Db = 0.0263 * exp(-28200/(RR*T0));      % Pawel
% Da = 2.36391687423300e-08;
% Db = 2.75733398155866e-07;
step_Da = 0.1*Da;
step_Db = 0.1*Db;
residual_idx = 8;

D_mat = zeros(5,2);
diff_LS_mat = zeros(5,7);

DaDb = zeros(max_iter, 2);
DaDb(1,:) = [Da Db];

diff_LS = zeros(max_iter, 10);

C4_oxidation;

load(strcat('output_file/', file_name));
C4_exp_data;
C4_least_square;

diff_LS_mat(1,1:3) = [wg_diff(2,1) d_diff(2,1) alpha_diff(2,1)];
sum_diff_LS = wg_diff(2,1) + d_diff(2,1) + alpha_diff(2,1);
weight_wg = wg_diff(2,1)/sum_diff_LS;
weight_d = d_diff(2,1)/sum_diff_LS;
weight_a = alpha_diff(2,1)/sum_diff_LS;
diff_LS_mat(1,4:5) = [diff_LS_mat(1,1)+diff_LS_mat(1,2)+diff_LS_mat(1,3) weight_wg*diff_LS_mat(1,1)+weight_d*diff_LS_mat(1,2)+weight_a*diff_LS_mat(1,3)];

if exist('wg_exp_UW', 'var') == 1
    diff_LS_mat(1,6:8) = [wg_diff_C4_CPUW_exp d_diff_C4_CPUW_exp alpha_diff_C4_CPUW_exp];
else
    diff_LS_mat(1,6:8) = [wg_diff_C4_CP_exp d_diff_C4_CP_exp alpha_diff_C4_CP_exp];
end
sum_diff_LS_exp = diff_LS_mat(1,6) + diff_LS_mat(1,7) + diff_LS_mat(1,8);
weight_wg_exp = diff_LS_mat(1,6)/sum_diff_LS_exp;
weight_d_exp = diff_LS_mat(1,7)/sum_diff_LS_exp;
weight_alpha_exp = diff_LS_mat(1,8)/sum_diff_LS_exp;
diff_LS_mat(1,9:10) = [0.2*diff_LS_mat(1,6)+0.2*diff_LS_mat(1,7)+0.6*diff_LS_mat(1,8) weight_wg_exp*diff_LS_mat(1,6)+weight_d_exp*diff_LS_mat(1,7)+weight_alpha_exp*diff_LS_mat(1,8)];

residual_sum = diff_LS_mat(1,residual_idx);

progress = waitbar(0, sprintf('step %d/%d', 0, max_iter), 'Name', 'Optimizing ...');

iter = 1;

while (iter <= max_iter) && (wg_diff_C4_CP > tolerance || d_diff_C4_CP > tolerance || alpha_diff_C4_CP > tolerance) && ((Da > Da_bound(1) && Da < Da_bound(2)) || (Db > Db_bound(1) && Db < Db_bound(2))) && (step_Da > 1e-3*Da && step_Db > 1e-3*Db)   %residual_sum > tolerance
    
    waitbar(iter / max_iter, progress, sprintf('step %d/%d', iter, max_iter));
    
    D_mat(1,:) = [Da Db];
    D_mat(2,:) = [min(Da+step_Da,Da_bound(2)) Db];
    D_mat(3,:) = [max(Da-step_Da,Da_bound(1)) Db];
    D_mat(4,:) = [Da min(Db+step_Db,Db_bound(2))];
    D_mat(5,:) = [Da max(Db-step_Db,Db_bound(1))];
    
    for iter_mat = 2:5
        Da = D_mat(iter_mat,1);
        Db = D_mat(iter_mat,2);
        
        C4_oxidation;
        
        load(strcat('output_file/', file_name));
        C4_exp_data;
        C4_least_square;
        
        diff_LS_mat(iter_mat,1:3) = [wg_diff(2,1) d_diff(2,1) alpha_diff(2,1)];
        sum_diff_LS = diff_LS_mat(iter_mat,1) + diff_LS_mat(iter_mat,2) + diff_LS_mat(iter_mat,3);
        weight_wg = diff_LS_mat(iter_mat,1)/sum_diff_LS;
        weight_d = diff_LS_mat(iter_mat,2)/sum_diff_LS;
        weight_a = diff_LS_mat(iter_mat,3)/sum_diff_LS;
        diff_LS_mat(iter_mat,4:5) = [sum_diff_LS weight_wg*diff_LS_mat(iter_mat,1)+weight_d*diff_LS_mat(iter_mat,2)+weight_a*diff_LS_mat(iter_mat,3)];
        
        if exist('wg_exp_UW', 'var') == 1
            diff_LS_mat(iter_mat,6:8) = [wg_diff_C4_CPUW_exp d_diff_C4_CPUW_exp alpha_diff_C4_CPUW_exp];
        else
            diff_LS_mat(iter_mat,6:8) = [wg_diff_C4_CP_exp d_diff_C4_CP_exp alpha_diff_C4_CP_exp];
        end
        sum_diff_LS_exp = diff_LS_mat(iter_mat,6) + diff_LS_mat(iter_mat,7) + diff_LS_mat(iter_mat,8);
        weight_wg_exp = diff_LS_mat(iter_mat,6)/sum_diff_LS_exp;
        weight_d_exp = diff_LS_mat(iter_mat,7)/sum_diff_LS_exp;
        weight_alpha_exp = diff_LS_mat(iter_mat,8)/sum_diff_LS_exp;
        diff_LS_mat(iter_mat,9:10) = [0.2*diff_LS_mat(iter_mat,6)+0.2*diff_LS_mat(iter_mat,7)+0.6*diff_LS_mat(iter_mat,8) weight_wg_exp*diff_LS_mat(iter_mat,6)+weight_d_exp*diff_LS_mat(iter_mat,7)+weight_alpha_exp*diff_LS_mat(iter_mat,8)];
    end
    
    [min_diff_LS_mat, index_min_diff_LS_mat] = min(diff_LS_mat(1:5,residual_idx));
    
    if min_diff_LS_mat < diff_LS_mat(1,residual_idx)
        Da = D_mat(index_min_diff_LS_mat,1);
        Db = D_mat(index_min_diff_LS_mat,2);
        diff_LS_mat(1,:) = diff_LS_mat(index_min_diff_LS_mat,:);
        DaDb(iter,:) = [Da Db];
        diff_LS(iter,:) = diff_LS_mat(index_min_diff_LS_mat,:);
    else
        step_Da = step_Da / 10;
        step_Db = step_Db / 10;
        Da = D_mat(1,1);
        Db = D_mat(1,2);
        DaDb(iter,:) = [Da Db];
        diff_LS(iter,:) = diff_LS_mat(1,:);
    end
    
    residual_sum = diff_LS(iter,residual_idx);
    iter = iter + 1;
    
    save(strcat('output_file/', file_name, '_HT_optimization_dicho'), 'D_mat', 'DaDb', 'Da_bound', 'Db_bound', 'step_Da', 'step_Db', '-regexp', 'diff_LS', 'diff_LS_mat', 'residual_idx', 'file_name');
    
end

close(progress);
clear progress;

toc


%% 2D-map in the Da - Db space
%{
C4_parameter;

N_Da = 21;
N_Db = 21;
N_D = N_Da * N_Db;

DaDb = zeros(N_D, 2);

diff_BJ_C4 = zeros(N_D, 2);
diff_Leist_C4 = zeros(N_D, 2);
diff_WPI_C4 = zeros(N_D, 3);
diff_Kawasaki_C4 = zeros(N_D, 3);
diff_CP_C4 = zeros(N_D, 3);

progress = waitbar(0, sprintf('step %d/%d', 0, N_D), 'Name', 'Optimizing ...');

for i_Da = 0:N_Da-1
    for i_Db = 0:N_Db-1
        
        i_DaDb = N_Da * i_Da + i_Db + 1;
        
        Da = 5e-9 + i_Da * 1e-8;
        Db = 5e-8 + i_Db * 1e-7;
        
        waitbar(i_DaDb / N_D, progress, sprintf('step %d/%d', i_DaDb, N_D));
        
        C4_oxidation;
        
        load(strcat('output_file/', file_name));
        C4_exp_data;
        C4_least_square;
        
        DaDb(i_DaDb, :) = [Da Db];

        diff_BJ_C4(i_DaDb, :) = [sqrt(wg_diff_BJ_C4) sqrt(d_diff_BJ_C4)];
        diff_Leist_C4(i_DaDb, :) = [sqrt(wg_diff_Leist_C4) sqrt(d_diff_Leist_C4)];
        diff_WPI_C4(i_DaDb, :) = [sqrt(wg_diff_WPI_C4) sqrt(d_diff_WPI_C4) sqrt(alpha_diff_WPI_C4)];
        diff_Kawasaki_C4(i_DaDb, :) = [sqrt(wg_diff_Kawasaki_C4) sqrt(d_diff_Kawasaki_C4) sqrt(alpha_diff_Kawasaki_C4)];
        diff_CP_C4(i_DaDb, :) = [sqrt(wg_diff_CP_C4) sqrt(d_diff_CP_C4) sqrt(alpha_diff_CP_C4)];
        
    end
end

close(progress);
clear progress;
save(strcat('output_file/', file_name, 'HT_optimization_2'), 'DaDb', '-regexp', '^(diff)*');
%}