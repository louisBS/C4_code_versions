%% C4 model
%%Generate the plots of the data

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

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
set(0, 'DefaultFigurePosition', [0 0.035 1 0.885]);
set(0, 'DefaultFigurePaperType', 'usletter');
set(0, 'DefaultFigurePaperOrientation', 'landscape');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat('output_file/', file_name, '_HT_optimization_dicho'));

[diff_LS_plot,index] = sort(diff_LS(:,residual_idx));

DaDb(index(1),:);

nb_pts = size(DaDb,1);

scatter3(DaDb(index(1:nb_pts),1), DaDb(index(1:nb_pts),2), diff_LS_plot(1:nb_pts), 1000, linspace(1,10,length(diff_LS(1:nb_pts,residual_idx))), '.');
xlabel('Da [cm^2/s]');
ylabel('Db [cm^2/s]');
zlabel('least square residual');
colormap jet;
colorbar;