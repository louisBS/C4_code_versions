%% C4 model
%%Generate the plots of the data

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated: 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set default graphic properties

set(0, 'DefaultFigureColor', 'White');
set(0, 'DefaultAxesFontSize', 20);      % 40 for presentation ; 30 else
set(0, 'DefaultAxesGridLineStyle', '-');
set(0, 'DefaultAxesXGrid', 'on');
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesXColor', 'Black');
set(0, 'DefaultAxesYColor', 'Black');
set(0, 'DefaultAxesZColor', 'Black');
set(0, 'DefaultAxeslineWidth', 1);
set(0, 'DefaultLineLineWidth', 4);      % 4 for presentation ; 3 else
set(0, 'DefaultLineMarkerSize', 15);    % 15 for presentation ; 10 else
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

%% Animated figure showing evolution of oxygen concentration with time
%%{

if strcmp(alloy, 'Zry4') && strcmp(temperature, '700-1200C')
    
    MovieObj = VideoWriter('animated_Co.avi');
    MovieObj.FrameRate = 10;
    open(MovieObj);
    
    fig_animated_Co = figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    x_plot = x * 1e4;
    x_ox_plot = x_ox * 1e4;
    xi_ox_plot = xi_ox * 1e4;
    xi_ab_plot = xi_ab * 1e4;
    d_protect_plot =  d_protect * 1e4;
    d_total_plot = d_total * 1e4;
    
    x_beta = [0 xi_ab_plot(1)];
    x_alpha = [xi_ab_plot(1) xi_ox_plot(1)];
    x_oxide_protect = [xi_ox_plot(1) xi_ox_plot(1)+d_protect_plot(1)];
    x_oxide_non_protect = [xi_ox_plot(1)+d_protect_plot(1) xi_ox_plot(1)+d_total_plot(1)];
    x_water = [xi_ox_plot(1)+d_total_plot(1) 1e4*e_Zr+d_total_plot(N+1)];
    y = [100 100];
    
    subplot(2,1,1);
    axis([500 1e4*e_Zr+d_total_plot(N+1) 0 1]);
    
    subplot(2,1,2);
    yyaxis right;
    hold on;
    alpha_phase = area([0 time_plot(N+1)], [850 850], 'FaceColor', [1 1 0], 'FaceAlpha', 0.3);
    beta_phase = fill([0 time_plot(N+1), time_plot(N+1) 0], [850 850, T(N+1)-Tk T(N+1)-Tk], 0, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.3);
    annotation('textbox', [0.75 0.15 0.12 0.05], 'String', '\alpha-Zr phase', 'Linewidth', 1.5, 'FontSize', 20, 'BackgroundColor', 'white');
    annotation('textbox', [0.75 0.23 0.12 0.12], 'String', '\beta-Zr phase + Oxygen-saturated \alpha-Zr phase', 'Linewidth', 1.5, 'FontSize', 20, 'BackgroundColor', 'white');
    yyaxis left;
    anim_d = animatedline('Color', 'r', 'LineWidth', 3);
    axis([0 time_plot(N+1) 0 d_total_plot(N+1)]);
    yyaxis right;
    anim_T = animatedline('Color', 'b', 'LineWidth', 3);
    
    
    axis([0 time_plot(N+1) T(N+1)-Tk T(1)-Tk]);
    
    
    for k = 1:N+1
        if mod(k, 10) == 1
            subplot(2,1,1);
            x_beta = [0 xi_ab_plot(k)];
            x_alpha = [xi_ab_plot(k) xi_ox_plot(k)];
            x_oxide_protect = [xi_ox_plot(k) xi_ox_plot(k)+d_protect_plot(k)];
            x_oxide_non_protect = [xi_ox_plot(k)+d_protect_plot(k) xi_ox_plot(k)+d_total_plot(k)];
            x_water = [xi_ox_plot(k)+d_total_plot(k) 1e4*e_Zr+d_total_plot(N+1)];
            
            if T(k) < T_transf
                alpha = area(x_beta, y, 'FaceColor', [1 1 0], 'FaceAlpha', 0.3);
                hold on;
                beta = area([0 0], [0 0], 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.3);
                alpha_O = area([0 0], [0 0], 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.6);
            else
                alpha = area([0 0], [0 0], 'FaceColor', [1 1 0], 'FaceAlpha', 0.3);
                hold on;
                beta = area(x_beta, y, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.3);
                alpha_O = area(x_alpha, y, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.6);
            end
            oxide_protect = area(x_oxide_protect, y, 'FaceColor', 'black', 'FaceAlpha', 0.3);
            oxide_non_protect = area(x_oxide_non_protect, y, 'FaceColor', 'black', 'FaceAlpha', 0.6);
            water = area(x_water, y, 'FaceColor', 'b', 'FaceAlpha', 0.3);
            
            O_content = plot(x_plot(:,k), 100*xo(:,k), 'Color', 'r', 'LineWidth', 3);
            plot([xi_ox_plot(k) xi_ox_plot(k)], 100*[xo_a_ox xo_ox(1,k)], 'r', 'LineWidth', 3);
            plot(xi_ox_plot(k)+x_ox_plot(:, k), 100*xo_ox(:, k), 'r', 'LineWidth', 3);
            hold off;
            
            axis([500 1e4*e_Zr+d_total_plot(N+1) 0 100]);
            title(sprintf('t = %.0f s', time_plot(k)));
            xlabel('cladding thickness [\mum]');
            ylabel('Oxygen concentration [at. %]');
            legend([O_content alpha beta alpha_O oxide_protect oxide_non_protect water], 'Oxygen concentration', '\alpha-Zr phase', '\beta-Zr phase', 'Oxygen-saturated \alpha-Zr phase', 'Protective oxide ZrO_2', 'Non protective oxide ZrO_2', 'water', 'Location', 'NorthWest');
            
            subplot(2,1,2);
            
            yyaxis left;
            addpoints(anim_d, time_plot(k), d_total_plot(k));
            xlabel(strcat('Exposure time [', time_unit, ']'));
            ylabel('Oxide thickness [\mum]');
            set(gca, 'YColor', 'r');
            set(gca, 'GridColor', 'k');
            
            yyaxis right;
            addpoints(anim_T, time_plot(k), T(k) - Tk);
            ylabel(sprintf('Temperature [%cC]', char(176)));
            set(gca, 'YColor', 'b');
            set(gca, 'GridColor', 'k');
            
            
            drawnow();
            
            frame = getframe(fig_animated_Co);
            writeVideo(MovieObj, frame);
        end
    end
    
    close(MovieObj);
    close(fig_animated_Co);
end
%}

%% Isotherm video

%{
if strcmp(alloy, 'Zry4')
    
    MovieObj = VideoWriter('animated_Co_1400C_100s.avi');
    MovieObj.FrameRate = 10;
    open(MovieObj);
    
    fig_animated_Co = figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    x_plot = x * 1e4;
    x_ox_plot = x_ox * 1e4;
    xi_ox_plot = xi_ox * 1e4;
    xi_ab_plot = xi_ab * 1e4;
    duct_b_plot = duct_b * 1e4;
    d_protect_plot =  d_protect * 1e4;
    d_total_plot = d_total * 1e4;
    
    x_duct_beta = [0 duct_b_plot(1)];
    x_brit_beta = [duct_b_plot(1) xi_ab_plot(1)];
    x_alpha = [xi_ab_plot(1) xi_ox_plot(1)];
    x_oxide_protect = [xi_ox_plot(1) xi_ox_plot(1)+d_protect_plot(1)];
    x_oxide_non_protect = [xi_ox_plot(1)+d_protect_plot(1) xi_ox_plot(1)+d_total_plot(1)];
    x_water = [xi_ox_plot(1)+d_total_plot(1) 1e4*e_Zr+d_total_plot(N+1)];
    y = [100 100];
    
    subplot(2,1,1);
    axis([300 1e4*e_Zr+d_total_plot(N+1) 0 1]);
    
    subplot(2,1,2);
    yyaxis left;
    anim_d = animatedline('Color', 'r', 'LineWidth', 3);
    axis([0 time_plot(N+1) 0 100]);
    yyaxis right;
    anim_T = animatedline('Color', 'b', 'LineWidth', 3);
    axis([0 time_plot(N+1) 0 100]);
    
    
    for k = 1:N+1
        if mod(k, 1) == 0
            subplot(2,1,1);
            x_duct_beta = [0 duct_b_plot(k)];
            x_brit_beta = [duct_b_plot(k) xi_ab_plot(k)];
            x_alpha = [xi_ab_plot(k) xi_ox_plot(k)];
            x_oxide_protect = [xi_ox_plot(k) xi_ox_plot(k)+d_protect_plot(k)];
            x_oxide_non_protect = [xi_ox_plot(k)+d_protect_plot(k) xi_ox_plot(k)+d_total_plot(k)];
            x_water = [xi_ox_plot(k)+d_total_plot(k) 1e4*e_Zr+d_total_plot(N+1)];
            
            duct_beta = area(x_duct_beta, y, 'FaceColor', [0 1 0], 'FaceAlpha', 0.6);
            hold on;
            brit_beta = area(x_brit_beta, y, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.6);
            alpha_O = area(x_alpha, y, 'FaceColor', [0 0.3 0], 'FaceAlpha', 0.6);

            oxide_protect = area(x_oxide_protect, y, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6);
            water = area(x_water, y, 'FaceColor', 'b', 'FaceAlpha', 0.6);
            
            O_content = plot(x_plot(:,k), 100*xo(:,k), 'Color', 'r', 'LineWidth', 3);
            plot([xi_ox_plot(k) xi_ox_plot(k)], 100*[xo_a_ox xo_ox(1,k)], 'r', 'LineWidth', 3);
            plot(xi_ox_plot(k)+x_ox_plot(:, k), 100*xo_ox(:, k), 'r', 'LineWidth', 3);
            hold off;
            
            axis([300 1e4*e_Zr+d_total_plot(N+1) 0 100]);
            title(sprintf('t = %.1f s', time_plot(k)));
            xlabel('cladding thickness [\mum]');
            ylabel('Oxygen concentration [at. %]');
            legend([O_content duct_beta brit_beta alpha_O oxide_protect water], 'Oxygen concentration', 'ductile \beta-Zr phase', 'brittle \beta-Zr phase', 'Oxygen-saturated \alpha-Zr phase', 'Protective oxide ZrO_2', 'water', 'Location', 'NorthWest');
            
            subplot(2,1,2);
            
            yyaxis left;
            addpoints(anim_d, time_plot(k), d_total_plot(k));
            axis([]);
            xlabel(strcat('Exposure time [', time_unit, ']'));
            ylabel('Oxide thickness [\mum]');
            set(gca, 'YColor', 'r');
            set(gca, 'GridColor', 'k');
            
            yyaxis right;
            addpoints(anim_T, time_plot(k), 100 * duct_b(k) / e_Zr);
            ylabel('Ductile \beta-Zr phase thickness [%]');
            set(gca, 'YColor', 'b');
            set(gca, 'GridColor', 'k');
            
            
            drawnow();
            
            frame = getframe(fig_animated_Co);
            writeVideo(MovieObj, frame);
        end
    end
    
    close(MovieObj);
    close(fig_animated_Co);
end
%}

%% Animated figure showing evolution of species concentration and the amount of oxidized Nb needed in the oxide with time
%{

MovieObj = VideoWriter('animated_C_ox.avi');
MovieObj.FrameRate=20;
open(MovieObj)

fig_animated_C_ox = figure();

subplot(2,1,1);
axis([0 100 1e17 1e23]);

subplot(2,1,2);
axis([0 100 0 3])


for k = 1:N+1
    subplot(2,1,1);
    
    semilogy(100*x_ox(:, k)/x_ox(end, k), Cv_ox(:, k), '-', 'Color', [1 0.5 0], 'lineWidth', 2);
    hold on
    semilogy(100*x_ox(:, k)/x_ox(end, k), Ce_ox(:, k), '-b', 'lineWidth', 2);
    hold on
    semilogy(100*x_ox(:, k)/x_ox(end, k), Ch_ox(:, k)', '-', 'Color', [0 0.6 0], 'lineWidth', 2);
    hold on
    semilogy(100*x_ox(:, k)/x_ox(end, k), Co_ox(:, k)', '-r', 'lineWidth', 2);
    hold off
    
    axis([0 100 1e17 1e23]);
    title(sprintf('t = %.1f days', time_plot(k)));
    xlabel('Normalized oxide thickness (%)', 'FontSize', 20);
    ylabel('Concentration (atom/cm^3)', 'FontSize', 20);
    legend('Vacancies concentration', 'Electron concentration', 'Hydrogen concentration', 'Oxygen concentration', 'Location', 'SouthWest');
    set(gca, 'FontSize', 20);
    grid on
    
    subplot(2,1,2);
    load(strcat('output_file/', file_name_C4), 'Nb_needed');
    plot(100*x_ox(:, t_plot)/x_ox(end, t_plot), 100*Nb_needed(:, t_plot), '-b', 'lineWidth', 1.5);
    hold on
    load(strcat('output_file/', file_name_C4_H), 'Nb_needed');
    plot(100*x_ox(:, t_plot)/x_ox(end, t_plot), 100*Nb_needed(:, t_plot), '-r', 'lineWidth', 1.5);
    hold on
    plot([0 100], [0.5 0.5], '--', 'Color', 'black', 'LineWidth', 1.5);
    hold off
    xlabel('Normalized oxide thickness (%)', 'FontSize', 20);
    ylabel('Oxidized Nb needed (wt %)', 'FontSize', 20);
    grid on
    
    drawnow();
    
    frame = getframe(fig_animated_C_ox);
    writeVideo(MovieObj, frame);
    
end

close(MovieObj);
close(fig_animated_C_ox);
%}

%% Animated figure showing the ratio of hydrogen and electron flux
%{

fig_animated_flux = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
figure(fig_animated_flux);

x_e = [0 d(1)];
y_e = [1 1];
x_h = [d(1) d(1)];
y_h = [1 1];

flux_e = area(x_e, y_e, 'FaceColor', 'b');
hold on
flux_h = area(x_h, y_h, 'FaceColor', [0 0.6 0]);
title('t = 0 days');
axis([0 0.001 0 1])

flux_movie = getframe(gcf);
[im, map] = rgb2ind(flux_movie.cdata, 256, 'nodither');
im(1,1,1,20) = 0;
frame_number = 1;
pause(0.01)

for k = 2:N+1
    x_e = [0 (1-fh_inst(k))*d(k)];
    x_h = [(1-fh_inst(k))*d(k) d(k)];
    set(flux_e, 'XData', x_e);
    set(flux_e, 'YData', y_e);
    set(flux_h, 'XData', x_h);
    set(flux_h, 'YData', y_h);
    title(sprintf('t = %d days', k-1));
    drawnow();
    
    flux_movie = getframe(gcf);
    frame_number = frame_number + 1;
    im(:, :, 1, frame_number) = rgb2ind(flux_movie.cdata, map, 'nodither');
    pause(0.01);
end

imwrite(im, map, 'animated_flux_ratio.gif', 'Delaytime_plot', 0, 'LoopCount', inf)
%}