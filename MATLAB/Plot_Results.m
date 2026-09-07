function figures = Plot_Results(results)
%PLOT_RESULTS Build and export all Task 1--6 figures.
% Task-specific plotters are local functions in this file.

    figures = gobjects(0);

    if nargin < 1 || ~isfield(results,'task1') || isempty(results.task1)
        return;
    end

    task1 = results.task1;
    if numel(task1) ~= 2
        error('Plot_Results:Task1CaseCount', ...
            'Task 1 requires step and multisine-command results.');
    end

    case_colors = lines(2);
    case_styles = {'-','--'};
    case_labels = {task1.label};
    time_label = '$t~[\mathrm{s}]$';

    % Tracking error
    state_titles = {'Roll-angle error','Pitch-angle error', ...
        'Yaw-angle error','Roll-rate error','Pitch-rate error', ...
        'Yaw-rate error'};
    y_labels = {'$e_{\phi}~[\mathrm{deg}]$', ...
        '$e_{\theta}~[\mathrm{deg}]$', ...
        '$e_{\psi}~[\mathrm{deg}]$', ...
        '$e_p~[\mathrm{rad/s}]$', ...
        '$e_q~[\mathrm{rad/s}]$', ...
        '$e_r~[\mathrm{rad/s}]$'};
    state_scale = [180/pi,180/pi,180/pi,1,1,1];

    figures(end+1) = figure('Color','w', ...
        'Name','Task 1 Tracking Error', ...
        'Position',[80 50 1080 740]);

    for state = 1:6
        subplot(3,2,state); hold on;
        for ic = 1:2
            plot(task1(ic).t,task1(ic).e(:,state)*state_scale(state), ...
                'Color',case_colors(ic,:), ...
                'LineStyle',case_styles{ic}, ...
                'LineWidth',1.55);
        end
        yline(0,'k:','LineWidth',0.9);
        latex_style(y_labels{state},time_label,state_titles{state});

        if state == 1
            legend(case_labels,'Interpreter','latex', ...
                'Location','best','FontSize',9);
        end
    end
    latex_super_title('Task 1: tracking error');

    % Control effort
    figures(end+1) = figure('Color','w', ...
        'Name','Task 1 Control Effort', ...
        'Position',[80 150 1200 390]);

    active_channels = [2,3,4];
    for input_index = 1:3
        subplot(1,3,input_index); hold on;
        for ic = 1:2
            plot(task1(ic).t,task1(ic).u(:,input_index), ...
                'Color',case_colors(ic,:), ...
                'LineStyle',case_styles{ic}, ...
                'LineWidth',1.55);
        end

        latex_style(sprintf('$u_%d$',active_channels(input_index)), ...
            time_label,sprintf('Active input $u_%d$', ...
            active_channels(input_index)));

        if input_index == 1
            legend(case_labels,'Interpreter','latex', ...
                'Location','best','FontSize',9);
        end
    end
    latex_super_title('Task 1: control effort');

    % Gain and uncertainty estimates
    figures(end+1) = matrix_estimate_figure(task1,case_colors, ...
        case_styles,case_labels,'Kx_hat','Kx_star',6,3, ...
        'Task 1 State Feedback Gain Estimates', ...
        'Task 1: state-feedback gain estimates','K_x');

    figures(end+1) = matrix_estimate_figure(task1,case_colors, ...
        case_styles,case_labels,'Kr_hat','Kr_star',3,3, ...
        'Task 1 Command Gain Estimates', ...
        'Task 1: command-gain estimates','K_r');

    figures(end+1) = matrix_estimate_figure(task1,case_colors, ...
        case_styles,case_labels,'Theta_hat','Theta_star',6,3, ...
        'Task 1 Uncertainty Estimates', ...
        'Task 1: uncertainty estimates','\Theta');

    % Parameter norms
    figures(end+1) = parameter_norm_figure(task1,case_colors, ...
        case_styles,case_labels);

    % Export PNG and EPS copies.
    if isfield(results,'output_dir') && ~isempty(results.output_dir)
        export_dir = fullfile(results.output_dir,'plots_eps');
        if ~exist(export_dir,'dir')
            mkdir(export_dir);
        end

        drawnow;
        previous_warning = warning('query', ...
            'MATLAB:print:ContentTypeImageSuggested');
        warning('off','MATLAB:print:ContentTypeImageSuggested');
        cleanup_warning = onCleanup(@() warning(previous_warning.state, ...
            'MATLAB:print:ContentTypeImageSuggested'));

        used_names = cell(1,numel(figures));
        for k = 1:numel(figures)
            fig = figures(k);
            base_name = regexprep(fig.Name,'[^a-zA-Z0-9]+','_');
            base_name = regexprep(base_name,'^_+|_+$','');
            if isempty(base_name)
                base_name = sprintf('Figure_%02d',fig.Number);
            end

            candidate_name = base_name;
            duplicate_index = 2;
            while any(strcmp(used_names,candidate_name))
                candidate_name = sprintf('%s_%02d', ...
                    base_name,duplicate_index);
                duplicate_index = duplicate_index+1;
            end
            used_names{k} = candidate_name;

            eps_path = fullfile(export_dir,[candidate_name '.eps']);
            png_path = fullfile(export_dir,[candidate_name '.png']);

            try
                exportgraphics(fig,eps_path, ...
                    'ContentType','vector','BackgroundColor','w');
                exportgraphics(fig,png_path, ...
                    'Resolution',300,'BackgroundColor','w');
            catch
                print(fig,fullfile(export_dir,candidate_name), ...
                    '-depsc','-r300');
                print(fig,png_path,'-dpng','-r300');
            end
        end
    end

    if isfield(results,'task2') && ~isempty(results.task2)
        task2_figures = Plot_Task2_Results(results);
        figures = [figures,task2_figures];
    end
    if isfield(results,'task3') && ~isempty(results.task3)
        task3_figures = Plot_Task3_Results(results);
        figures = [figures,task3_figures];
    end
    if isfield(results,'task4') && ~isempty(results.task4)
        task4_figures = Plot_Task4_Results(results);
        figures = [figures,task4_figures];
    end
    if isfield(results,'task5') && ~isempty(results.task5)
        task5_figures = Plot_Task5_Results(results);
        figures = [figures,task5_figures];
    end
    if isfield(results,'task6') && ~isempty(results.task6)
        task6_figures = Plot_Task6_Results(results);
        figures = [figures,task6_figures];
    end
end

function fig = matrix_estimate_figure(task1,case_colors,case_styles, ...
        case_labels,history_field,ideal_field,n_rows,n_cols, ...
        figure_name,super_title,symbol)
% Plot every estimated matrix entry against its ideal value.

    figure_height = max(650,175*n_rows);
    fig = figure('Color','w','Name',figure_name, ...
        'Position',[40 30 1250 figure_height]);

    for row = 1:n_rows
        for column = 1:n_cols
            tile = (row-1)*n_cols+column;
            subplot(n_rows,n_cols,tile); hold on;

            plot_handles = gobjects(1,2);
            for ic = 1:2
                history = task1(ic).(history_field);
                trajectory = squeeze(history(:,row,column));
                plot_handles(ic) = plot(task1(ic).t,trajectory, ...
                    'Color',case_colors(ic,:), ...
                    'LineStyle',case_styles{ic}, ...
                    'LineWidth',1.25);
            end

            ideal_matrix = task1(1).(ideal_field);
            ideal_handle = yline(ideal_matrix(row,column),'k:', ...
                'LineWidth',1.1);

            % sprintf interprets "\w" in "\widehat" as an invalid escape.
            title_text = ['$\widehat{' symbol '}_{' ...
                num2str(row) num2str(column) '}$'];

            x_label = '';
            if row == n_rows
                x_label = '$t~[\mathrm{s}]$';
            end

            y_label = '';
            if column == 1
                y_label = 'Estimate';
            end

            latex_style(y_label,x_label,title_text);

            if row == 1 && column == 1
                legend([plot_handles,ideal_handle], ...
                    [case_labels,{'Ideal value'}], ...
                    'Interpreter','latex','Location','best', ...
                    'FontSize',7.5);
            end
        end
    end

    latex_super_title(super_title);
end

function fig = parameter_norm_figure(task1,case_colors,case_styles,case_labels)
% Summarize estimate sizes and estimation errors.

    history_fields = {'Kx_hat','Kr_hat','Theta_hat'};
    ideal_fields = {'Kx_star','Kr_star','Theta_star'};
    error_fields = {'Kx_error_norm','Kr_error_norm','Theta_error_norm'};
    symbols = {'K_x','K_r','\Theta'};
    panel_titles = {'State-feedback gain','Command gain', ...
        'Uncertainty parameter'};

    fig = figure('Color','w', ...
        'Name','Task 1 Adaptive Parameter Estimates', ...
        'Position',[80 50 1200 680]);

    for parameter_index = 1:3
        subplot(2,3,parameter_index); hold on;
        plot_handles = gobjects(1,2);
        for ic = 1:2
            history = task1(ic).(history_fields{parameter_index});
            estimate_norm = sqrt(sum(reshape(history,size(history,1),[]).^2,2));
            plot_handles(ic) = plot(task1(ic).t,estimate_norm, ...
                'Color',case_colors(ic,:), ...
                'LineStyle',case_styles{ic}, ...
                'LineWidth',1.55);
        end

        ideal_norm = norm(task1(1).(ideal_fields{parameter_index}),'fro');
        ideal_handle = yline(ideal_norm,'k:','LineWidth',1.1);
        latex_style( ...
            ['$\|\widehat{' symbols{parameter_index} '}(t)\|_F$'], ...
            '',panel_titles{parameter_index});

        if parameter_index == 1
            legend([plot_handles,ideal_handle], ...
                [case_labels,{'Ideal norm'}], ...
                'Interpreter','latex','Location','best','FontSize',8);
        end

        subplot(2,3,parameter_index+3); hold on;
        for ic = 1:2
            plot(task1(ic).t,task1(ic).(error_fields{parameter_index}), ...
                'Color',case_colors(ic,:), ...
                'LineStyle',case_styles{ic}, ...
                'LineWidth',1.55);
        end
        latex_style( ...
            ['$\|\widetilde{' symbols{parameter_index} '}(t)\|_F$'], ...
            '$t~[\mathrm{s}]$',[panel_titles{parameter_index} ' error']);
    end

    latex_super_title('Task 1: adaptive parameter norms and errors');
end

function latex_style(y_label,x_label,title_text)
% Apply the common LaTeX plot style.

    grid on;
    grid minor;
    box on;
    set(gca,'FontSize',10,'LineWidth',1.0, ...
        'TickLabelInterpreter','latex');

    if ~isempty(y_label)
        ylabel(y_label,'Interpreter','latex','FontSize',11);
    end
    if ~isempty(x_label)
        xlabel(x_label,'Interpreter','latex','FontSize',11);
    end
    if nargin >= 3 && ~isempty(title_text)
        title(title_text,'Interpreter','latex','FontSize',11, ...
            'FontWeight','bold');
    end
end

function latex_super_title(title_text)
% Apply a consistent figure-level title.

    if exist('sgtitle','file') == 2
        sgtitle(title_text,'Interpreter','latex','FontSize',15, ...
            'FontWeight','bold');
    else
        annotation('textbox',[0,0.955,1,0.04], ...
            'String',title_text,'Interpreter','latex', ...
            'EdgeColor','none','HorizontalAlignment','center', ...
            'FontSize',15,'FontWeight','bold');
    end
end


function figures = Plot_Task2_Results(results)
% Build the Task 2 comparison figures.

    figures = gobjects(0);
    if nargin < 1 || ~isfield(results,'task2')
        return;
    end

    direct = results.task2.direct;
    indirect = results.task2.indirect;
    if numel(direct) ~= 2 || numel(indirect) ~= 2
        error('Plot_Task2_Results:CaseCount', ...
            'Task 2 requires step and multisine cases for both controllers.');
    end

    controller_colors = lines(2);
    controller_labels = {'Direct MRAC','Indirect MRAC'};
    time_label = '$t~[\mathrm{s}]$';

    % Tracking and control effort
    figures(end+1) = figure('Color','w', ...
        'Name','Task 2 Transient Comparison', ...
        'Position',[70 50 1200 720]);

    for command_index = 1:2
        subplot(2,2,command_index); hold on;
        plot(direct(command_index).t, ...
            sqrt(sum(direct(command_index).e.^2,2)), ...
            'Color',controller_colors(1,:),'LineWidth',1.55);
        plot(indirect(command_index).t, ...
            sqrt(sum(indirect(command_index).e.^2,2)), ...
            '--','Color',controller_colors(2,:),'LineWidth',1.55);
        latex_style('$\|e(t)\|_2$','', ...
            [direct(command_index).label ': tracking']);
        if command_index == 1
            legend(controller_labels,'Interpreter','latex', ...
                'Location','best','FontSize',9);
        end

        subplot(2,2,command_index+2); hold on;
        plot(direct(command_index).t, ...
            sqrt(sum(direct(command_index).u.^2,2)), ...
            'Color',controller_colors(1,:),'LineWidth',1.55);
        plot(indirect(command_index).t, ...
            sqrt(sum(indirect(command_index).u.^2,2)), ...
            '--','Color',controller_colors(2,:),'LineWidth',1.55);
        latex_style('$\|u(t)\|_2$',time_label, ...
            [direct(command_index).label ': control effort']);
    end
    latex_super_title('Task 2: direct versus indirect MRAC transients');

    % Parameter and certainty-equivalence errors
    figures(end+1) = figure('Color','w', ...
        'Name','Task 2 Parameter Convergence', ...
        'Position',[50 30 1300 760]);

    for command_index = 1:2
        row_offset = 3*(command_index-1);

        subplot(2,3,row_offset+1); hold on;
        plot(indirect(command_index).t, ...
            indirect(command_index).Lambda_error_norm, ...
            'Color',controller_colors(2,:),'LineWidth',1.55);
        latex_style('$\|\widehat{\Lambda}-\Lambda\|_F$', ...
            conditional_time_label(command_index,time_label), ...
            [indirect(command_index).label ': effectiveness']);
        include_zero_lower_bound();
        if command_index == 1
            legend({'Indirect estimate'},'Interpreter','latex', ...
                'Location','best','FontSize',8);
        end

        subplot(2,3,row_offset+2); hold on;
        plot(direct(command_index).t, ...
            direct(command_index).Theta_error_norm, ...
            'Color',controller_colors(1,:),'LineWidth',1.55);
        plot(indirect(command_index).t, ...
            indirect(command_index).Theta_error_norm, ...
            '--','Color',controller_colors(2,:),'LineWidth',1.55);
        latex_style('$\|\widetilde{\Theta}\|_F$', ...
            conditional_time_label(command_index,time_label), ...
            [indirect(command_index).label ': uncertainty']);
        include_zero_lower_bound();
        if command_index == 1
            legend(controller_labels,'Interpreter','latex', ...
                'Location','best','FontSize',8);
        end

        subplot(2,3,row_offset+3); hold on;
        direct_gain_error = sqrt( ...
            direct(command_index).Kx_error_norm.^2 ...
            + direct(command_index).Kr_error_norm.^2);
        indirect_gain_error = sqrt( ...
            indirect(command_index).Kx_error_norm.^2 ...
            + indirect(command_index).Kr_error_norm.^2);
        plot(direct(command_index).t,direct_gain_error, ...
            'Color',controller_colors(1,:),'LineWidth',1.55);
        plot(indirect(command_index).t,indirect_gain_error, ...
            '--','Color',controller_colors(2,:),'LineWidth',1.55);
        latex_style( ...
            '$\sqrt{\|\widetilde{K}_x\|_F^2+\|\widetilde{K}_r\|_F^2}$', ...
            conditional_time_label(command_index,time_label), ...
            [indirect(command_index).label ': matched gains']);
        include_zero_lower_bound();
    end
    latex_super_title('Task 2: parameter and gain convergence');

    % Predictor error and finite-window excitation
    figures(end+1) = figure('Color','w', ...
        'Name','Task 2 Identification and PE', ...
        'Position',[70 50 1200 720]);

    for command_index = 1:2
        subplot(2,2,command_index); hold on;
        plot(indirect(command_index).t, ...
            indirect(command_index).prediction_error_norm, ...
            'Color',controller_colors(2,:),'LineWidth',1.55);
        latex_style('$\|x-\widehat{x}_{\mathrm{id}}\|_2$','', ...
            [indirect(command_index).label ': predictor error']);

        subplot(2,2,command_index+2); hold on;
        semilogy(indirect(command_index).t, ...
            max(indirect(command_index).pe_min_eigenvalue,eps), ...
            'Color',controller_colors(2,:),'LineWidth',1.55);
        set(gca,'YScale','log');
        latex_style( ...
            '$\lambda_{\min}$ of windowed regressor Gramian', ...
            time_label,[indirect(command_index).label ': PE indicator']);
    end
    latex_super_title( ...
        'Task 2: identifier prediction error and excitation');

    if isfield(results,'output_dir') && ~isempty(results.output_dir)
        export_task2_figures(figures,results.output_dir);
    end
end

function include_zero_lower_bound()
% Include zero so norm convergence is not visually exaggerated.

    current_limits = ylim;
    ylim([0,current_limits(2)]);
end

function label = conditional_time_label(command_index,time_label)
    if command_index == 2
        label = time_label;
    else
        label = '';
    end
end

function export_task2_figures(figures,output_dir)
    export_dir = fullfile(output_dir,'plots_eps');
    if ~exist(export_dir,'dir')
        mkdir(export_dir);
    end

    drawnow;
    previous_warning = warning('query', ...
        'MATLAB:print:ContentTypeImageSuggested');
    warning('off','MATLAB:print:ContentTypeImageSuggested');
    cleanup_warning = onCleanup(@() warning(previous_warning.state, ...
        'MATLAB:print:ContentTypeImageSuggested'));

    used_names = cell(1,numel(figures));
    for k = 1:numel(figures)
        fig = figures(k);
        base_name = regexprep(fig.Name,'[^a-zA-Z0-9]+','_');
        base_name = regexprep(base_name,'^_+|_+$','');
        if isempty(base_name)
            base_name = sprintf('Task_2_Figure_%02d',fig.Number);
        end

        candidate_name = base_name;
        duplicate_index = 2;
        while any(strcmp(used_names,candidate_name))
            candidate_name = sprintf('%s_%02d', ...
                base_name,duplicate_index);
            duplicate_index = duplicate_index+1;
        end
        used_names{k} = candidate_name;

        eps_path = fullfile(export_dir,[candidate_name '.eps']);
        png_path = fullfile(export_dir,[candidate_name '.png']);
        try
            exportgraphics(fig,eps_path, ...
                'ContentType','vector','BackgroundColor','w');
            exportgraphics(fig,png_path, ...
                'Resolution',300,'BackgroundColor','w');
        catch
            print(fig,fullfile(export_dir,candidate_name), ...
                '-depsc','-r300');
            print(fig,png_path,'-dpng','-r300');
        end
    end
end

function figures = Plot_Task3_Results(results)
% Build the Task 3 disturbance-rejection figures.

    figures = gobjects(0);
    if nargin < 1 || ~isfield(results,'task3')
        return;
    end

    baseline = results.task3.baseline;
    integral = results.task3.integral;
    controller_colors = lines(2);
    controller_labels = {baseline.controller,integral.controller};
    disturbance_time = integral.disturbance_start_time;
    time_label = '$t~[\mathrm{s}]$';

    % Output error, control, and integral state
    figures(end+1) = figure('Color','w', ...
        'Name','Task 3 Disturbance Rejection', ...
        'Position',[60 30 1250 820]);

    output_titles = {'Roll-output error','Pitch-output error', ...
        'Yaw-output error'};
    output_labels = {'$e_{\phi}~[\mathrm{deg}]$', ...
        '$e_{\theta}~[\mathrm{deg}]$', ...
        '$e_{\psi}~[\mathrm{deg}]$'};

    for output_index = 1:3
        subplot(3,2,2*output_index-1); hold on;
        plot(baseline.t,baseline.output_error(:,output_index)*180/pi, ...
            'Color',controller_colors(1,:),'LineWidth',1.55);
        plot(integral.t,integral.output_error(:,output_index)*180/pi, ...
            '--','Color',controller_colors(2,:),'LineWidth',1.55);
        disturbance_marker(disturbance_time);
        latex_style(output_labels{output_index},time_label, ...
            output_titles{output_index});
        if output_index == 1
            legend(controller_labels,'Interpreter','latex', ...
                'Location','best','FontSize',8);
        end
    end

    subplot(3,2,2); hold on;
    plot(baseline.t,vecnorm(baseline.output_error,2,2)*180/pi, ...
        'Color',controller_colors(1,:),'LineWidth',1.55);
    plot(integral.t,vecnorm(integral.output_error,2,2)*180/pi, ...
        '--','Color',controller_colors(2,:),'LineWidth',1.55);
    disturbance_marker(disturbance_time);
    latex_style('$\|e_y(t)\|_2~[\mathrm{deg}]$','', ...
        'Output-error norm');

    subplot(3,2,4); hold on;
    plot(baseline.t,vecnorm(baseline.u,2,2), ...
        'Color',controller_colors(1,:),'LineWidth',1.55);
    plot(integral.t,vecnorm(integral.u,2,2), ...
        '--','Color',controller_colors(2,:),'LineWidth',1.55);
    disturbance_marker(disturbance_time);
    latex_style('$\|u(t)\|_2$','', ...
        'Control effort');

    subplot(3,2,6); hold on;
    plot(integral.t,vecnorm(integral.eI,2,2)*180/pi, ...
        'Color',controller_colors(2,:),'LineWidth',1.55);
    disturbance_marker(disturbance_time);
    latex_style('$\|e_I(t)\|_2~[\mathrm{deg\,s}]$',time_label, ...
        'Integral-error state');
    latex_super_title( ...
        'Task 3: rejection of a matched constant input disturbance');

    % Integral-gain estimates
    figures(end+1) = figure('Color','w', ...
        'Name','Task 3 Integral Gain Estimates', ...
        'Position',[90 45 1100 760]);

    for row = 1:3
        for column = 1:3
            subplot(3,3,(row-1)*3+column); hold on;
            estimate = squeeze(integral.Ki_hat(:,row,column));
            plot(integral.t,estimate, ...
                'Color',controller_colors(2,:),'LineWidth',1.45);
            yline(integral.Ki_star(row,column),'k:', ...
                'LineWidth',1.1);
            disturbance_marker(disturbance_time);
            title_text = ['$\widehat{K}_{I,' ...
                num2str(row) num2str(column) '}$'];
            x_label = '';
            if row == 3
                x_label = time_label;
            end
            y_label = '';
            if column == 1
                y_label = 'Estimate';
            end
            latex_style(y_label,x_label,title_text);
            if row == 1 && column == 1
                legend({'Estimate','Ideal value'}, ...
                    'Interpreter','latex','Location','best','FontSize',8);
            end
        end
    end
    latex_super_title('Task 3: adaptive integral-gain estimates');

    if isfield(results,'output_dir') && ~isempty(results.output_dir)
        export_task3_figures(figures,results.output_dir);
    end
end

function disturbance_marker(disturbance_time)
    xline(disturbance_time,'k--','LineWidth',1.0);
end

function export_task3_figures(figures,output_dir)
    export_dir = fullfile(output_dir,'plots_eps');
    if ~exist(export_dir,'dir')
        mkdir(export_dir);
    end

    drawnow;
    previous_warning = warning('query', ...
        'MATLAB:print:ContentTypeImageSuggested');
    warning('off','MATLAB:print:ContentTypeImageSuggested');
    cleanup_warning = onCleanup(@() warning(previous_warning.state, ...
        'MATLAB:print:ContentTypeImageSuggested'));

    used_names = cell(1,numel(figures));
    for k = 1:numel(figures)
        fig = figures(k);
        base_name = regexprep(fig.Name,'[^a-zA-Z0-9]+','_');
        base_name = regexprep(base_name,'^_+|_+$','');
        if isempty(base_name)
            base_name = sprintf('Task_3_Figure_%02d',fig.Number);
        end

        candidate_name = base_name;
        duplicate_index = 2;
        while any(strcmp(used_names,candidate_name))
            candidate_name = sprintf('%s_%02d', ...
                base_name,duplicate_index);
            duplicate_index = duplicate_index+1;
        end
        used_names{k} = candidate_name;

        eps_path = fullfile(export_dir,[candidate_name '.eps']);
        png_path = fullfile(export_dir,[candidate_name '.png']);
        try
            exportgraphics(fig,eps_path, ...
                'ContentType','vector','BackgroundColor','w');
            exportgraphics(fig,png_path, ...
                'Resolution',300,'BackgroundColor','w');
        catch
            print(fig,fullfile(export_dir,candidate_name), ...
                '-depsc','-r300');
            print(fig,png_path,'-dpng','-r300');
        end
    end
end

function figures = Plot_Task4_Results(results)
% Compare the baseline with all four robust modifications.

    figures = gobjects(0);
    if nargin < 1 || ~isfield(results,'task4')
        return;
    end

    cases = results.task4.cases;
    n_cases = numel(cases);
    colors = lines(n_cases);
    line_styles = {'-','--','-.',':','-'};
    labels = cellfun(@(item) item.controller,cases, ...
        'UniformOutput',false);
    method_names = cellfun(@(item) item.robust_method,cases, ...
        'UniformOutput',false);
    labels{strcmp(method_names,'sigma')} = '$\sigma$-modification';
    labels{strcmp(method_names,'emod')} = '$e$-modification';
    time_label = '$t~[\mathrm{s}]$';

    % Main comparison
    figures(end+1) = figure('Color','w', ...
        'Name','Task 4 Robust Modification Comparison', ...
        'Position',[50 25 1260 820]);

    subplot(2,2,1); hold on;
    for case_index = 1:n_cases
        item = cases{case_index};
        plot(item.t,vecnorm(item.e,2,2), ...
            'Color',colors(case_index,:), ...
            'LineStyle',line_styles{case_index}, ...
            'LineWidth',1.45);
    end
    latex_style('$\|e(t)\|_2$',time_label,'Tracking-error norm');
    legend(labels,'Interpreter','latex','Location','best', ...
        'FontSize',8);

    subplot(2,2,2); hold on;
    for case_index = 1:n_cases
        item = cases{case_index};
        plot(item.t,vecnorm(item.u,2,2), ...
            'Color',colors(case_index,:), ...
            'LineStyle',line_styles{case_index}, ...
            'LineWidth',1.45);
    end
    latex_style('$\|u(t)\|_2$',time_label, ...
        'Control effort (five curves nearly overlap)');

    subplot(2,2,3); hold on;
    for case_index = 1:n_cases
        item = cases{case_index};
        plot(item.t,item.parameter_error_norm, ...
            'Color',colors(case_index,:), ...
            'LineStyle',line_styles{case_index}, ...
            'LineWidth',1.45);
    end
    latex_style('$\|\widetilde{\vartheta}(t)\|_2$', ...
        time_label,'Combined parameter-estimation error');

    subplot(2,2,4); hold on;
    for case_index = 1:n_cases
        item = cases{case_index};
        plot(item.t,item.parameter_estimate_norm, ...
            'Color',colors(case_index,:), ...
            'LineStyle',line_styles{case_index}, ...
            'LineWidth',1.45);
    end
    latex_style('$\|\widehat{\vartheta}(t)\|_2$', ...
        time_label,'Combined estimate norm');
    latex_super_title( ...
        'Task 4: robust MRAC under a common bounded disturbance');

    % Estimation errors and modification activity
    figures(end+1) = figure('Color','w', ...
        'Name','Task 4 Parameter Drift and Activity', ...
        'Position',[80 35 1260 820]);

    parameter_fields = {'Kx_error_norm','Kr_error_norm', ...
        'Theta_error_norm'};
    parameter_titles = {'State-feedback gain error', ...
        'Command-gain error','Uncertainty-parameter error'};
    parameter_symbols = {'K_x','K_r','\Theta'};

    for parameter_index = 1:3
        subplot(2,2,parameter_index); hold on;
        for case_index = 1:n_cases
            item = cases{case_index};
            plot(item.t,item.(parameter_fields{parameter_index}), ...
                'Color',colors(case_index,:), ...
                'LineStyle',line_styles{case_index}, ...
                'LineWidth',1.45);
        end
        latex_style( ...
            ['$\|\widetilde{' parameter_symbols{parameter_index} ...
             '}(t)\|_F$'],time_label, ...
            parameter_titles{parameter_index});
        if parameter_index == 1
            legend(labels,'Interpreter','latex','Location','best', ...
                'FontSize',8);
        end
    end

    subplot(2,2,4); hold on;
    deadzone_case = cases{strcmp(method_names,'deadzone')};
    projection_case = cases{strcmp(method_names,'projection')};
    stairs(deadzone_case.t,double(deadzone_case.adaptation_active), ...
        'Color',colors(2,:),'LineWidth',1.35);
    stairs(projection_case.t, ...
        double(projection_case.projection_active_count > 0), ...
        'Color',colors(5,:),'LineWidth',1.35);
    ylim([-0.05,1.05]);
    yticks([0,1]);
    yticklabels({'inactive','active'});
    latex_style('Modification activity',time_label, ...
        'Dead-zone adaptation and projection activity');
    legend({'Dead-zone adaptation','Projection boundary'}, ...
        'Interpreter','latex','Location','best','FontSize',8);
    latex_super_title('Task 4: parameter drift and robust-law activity');

    if isfield(results,'output_dir') && ~isempty(results.output_dir)
        export_task4_figures(figures,results.output_dir);
    end
end

function export_task4_figures(figures,output_dir)
    export_dir = fullfile(output_dir,'plots_eps');
    if ~exist(export_dir,'dir')
        mkdir(export_dir);
    end

    drawnow;
    previous_warning = warning('query', ...
        'MATLAB:print:ContentTypeImageSuggested');
    warning('off','MATLAB:print:ContentTypeImageSuggested');
    cleanup_warning = onCleanup(@() warning(previous_warning.state, ...
        'MATLAB:print:ContentTypeImageSuggested'));

    for figure_index = 1:numel(figures)
        fig = figures(figure_index);
        base_name = regexprep(fig.Name,'[^a-zA-Z0-9]+','_');
        base_name = regexprep(base_name,'^_+|_+$','');
        eps_path = fullfile(export_dir,[base_name '.eps']);
        png_path = fullfile(export_dir,[base_name '.png']);
        try
            exportgraphics(fig,eps_path, ...
                'ContentType','vector','BackgroundColor','w');
            exportgraphics(fig,png_path, ...
                'Resolution',300,'BackgroundColor','w');
        catch
            print(fig,fullfile(export_dir,base_name), ...
                '-depsc','-r300');
            print(fig,png_path,'-dpng','-r300');
        end
    end
end

function figures = Plot_Task5_Results(results)
% Build the LTR and output-feedback comparisons.

    figures = gobjects(0);
    if nargin < 1 || ~isfield(results,'task5')
        return;
    end

    design = results.task5.design;
    cases = results.task5.cases;
    n_rho = numel(design.rho);
    recovery_colors = lines(n_rho);
    case_colors = lines(numel(cases));
    case_styles = {'-','--','-.',':'};
    case_labels = cellfun(@(item) item.controller,cases, ...
        'UniformOutput',false);
    output_case_labels = case_labels(2:end);
    time_label = '$t~[\mathrm{s}]$';

    % Frequency-domain loop-transfer recovery
    figures(end+1) = figure('Color','w', ...
        'Name','Task 5 Loop Transfer Recovery', ...
        'Position',[45 25 1260 820]);

    singular_value_columns = [1,size(design.full_state_loop_sv,2)];
    singular_value_titles = { ...
        'Maximum loop singular value', ...
        'Minimum loop singular value'};
    for panel_index = 1:2
        subplot(2,2,panel_index); hold on;
        column = singular_value_columns(panel_index);
        loglog(design.frequency, ...
            design.full_state_loop_sv(:,column), ...
            'k','LineWidth',2.0);
        for rho_index = 1:n_rho
            loglog(design.frequency, ...
                design.observer_loop_sv(:,column,rho_index), ...
                'Color',recovery_colors(rho_index,:), ...
                'LineWidth',1.35);
        end
        set(gca,'XScale','log','YScale','log');
        latex_style('Singular value','$\omega~[\mathrm{rad/s}]$', ...
            singular_value_titles{panel_index});
        if panel_index == 1
            legend(recovery_legend(design.rho), ...
                'Interpreter','latex','Location','best', ...
                'FontSize',8);
        end
    end

    subplot(2,2,3);
    loglog(design.rho,design.recovery_log_sv_rmse, ...
        'o-','Color',[0.20,0.45,0.75], ...
        'MarkerFaceColor',[0.20,0.45,0.75], ...
        'LineWidth',1.55);
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    latex_style('Log-singular-value RMSE','$\rho$', ...
        'Recovery error (smaller is better)');

    subplot(2,2,4);
    loglog(design.rho,design.observer_gain_norm, ...
        's-','Color',[0.85,0.33,0.10], ...
        'MarkerFaceColor',[0.85,0.33,0.10], ...
        'LineWidth',1.55);
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    latex_style('$\|L(\rho)\|_2$','$\rho$', ...
        'Observer-gain cost of recovery');
    latex_super_title( ...
        'Task 5: observer-based loop-transfer recovery');

    % Closed-loop comparison
    figures(end+1) = figure('Color','w', ...
        'Name','Task 5 Output Feedback Comparison', ...
        'Position',[65 30 1260 820]);

    subplot(2,2,1); hold on;
    for case_index = 1:numel(cases)
        item = cases{case_index};
        plot(item.t,vecnorm(item.e,2,2), ...
            'Color',case_colors(case_index,:), ...
            'LineStyle',case_styles{case_index}, ...
            'LineWidth',1.45);
    end
    latex_style('$\|x-x_m\|_2$',time_label, ...
        'True tracking-error norm');
    legend(case_labels,'Interpreter','latex','Location','best', ...
        'FontSize',8);

    subplot(2,2,2); hold on;
    for case_index = 2:numel(cases)
        item = cases{case_index};
        plot(item.t,vecnorm(item.observer_error,2,2), ...
            'Color',case_colors(case_index,:), ...
            'LineStyle',case_styles{case_index}, ...
            'LineWidth',1.45);
    end
    set(gca,'YScale','log');
    latex_style('$\|x-\widehat{x}\|_2$',time_label, ...
        'Observer-error norm');
    legend(output_case_labels,'Interpreter','latex', ...
        'Location','best','FontSize',8);

    subplot(2,2,3); hold on;
    for case_index = 1:numel(cases)
        item = cases{case_index};
        plot(item.t,vecnorm(item.u,2,2), ...
            'Color',case_colors(case_index,:), ...
            'LineStyle',case_styles{case_index}, ...
            'LineWidth',1.45);
    end
    latex_style('$\|u(t)\|_2$',time_label,'Control effort');

    subplot(2,2,4); hold on;
    for case_index = 1:numel(cases)
        item = cases{case_index};
        plot(item.t,item.parameter_error_norm, ...
            'Color',case_colors(case_index,:), ...
            'LineStyle',case_styles{case_index}, ...
            'LineWidth',1.45);
    end
    latex_style('$\|\widetilde{\vartheta}(t)\|_2$', ...
        time_label,'Baseline MRAC parameter error');
    latex_super_title( ...
        'Task 5: full-state and output-feedback MRAC comparison');

    % Observer and augmentation diagnostics
    figures(end+1) = figure('Color','w', ...
        'Name','Task 5 Observer and Augmentation Diagnostics', ...
        'Position',[80 35 1280 800]);

    rate_names = {'Roll-rate','Pitch-rate','Yaw-rate'};
    rate_symbols = {'p','q','r'};
    for rate_index = 1:3
        subplot(2,3,rate_index); hold on;
        state_index = rate_index+3;
        for case_index = 2:numel(cases)
            item = cases{case_index};
            plot(item.t,item.observer_error(:,state_index)*180/pi, ...
                'Color',case_colors(case_index,:), ...
                'LineStyle',case_styles{case_index}, ...
                'LineWidth',1.35);
        end
        latex_style( ...
            ['$' rate_symbols{rate_index} ...
             '-\widehat{' rate_symbols{rate_index} ...
             '}~[\mathrm{deg/s}]$'],time_label, ...
            [rate_names{rate_index} ' estimation error']);
        if rate_index == 1
            legend(output_case_labels,'Interpreter','latex', ...
                'Location','best','FontSize',8);
        end
    end

    subplot(2,3,4); hold on;
    for case_index = 2:numel(cases)
        item = cases{case_index};
        plot(item.t,vecnorm(item.innovation,2,2)*180/pi, ...
            'Color',case_colors(case_index,:), ...
            'LineStyle',case_styles{case_index}, ...
            'LineWidth',1.35);
    end
    latex_style('$\|y-C\widehat{x}\|_2~[\mathrm{deg}]$', ...
        time_label,'Measured-output innovation');

    augmented_case = cases{strcmp(cellfun(@(item) ...
        item.configuration,cases,'UniformOutput',false),'augmented')};
    subplot(2,3,5);
    plot(augmented_case.t,augmented_case.W_error_norm, ...
        'Color',case_colors(4,:),'LineWidth',1.45);
    latex_style('$\|\widetilde{W}(t)\|_F$',time_label, ...
        'Augmentation parameter error');

    subplot(2,3,6);
    plot(augmented_case.t, ...
        vecnorm(augmented_case.u_augmentation,2,2), ...
        'Color',case_colors(4,:),'LineWidth',1.45);
    latex_style('$\|u_{\mathrm{aug}}(t)\|_2$',time_label, ...
        'Adaptive augmentation effort');
    latex_super_title( ...
        'Task 5: observer and adaptive-augmentation diagnostics');

    if isfield(results,'output_dir') && ~isempty(results.output_dir)
        export_task5_figures(figures,results.output_dir);
    end
end

function labels = recovery_legend(rho_values)
    labels = cell(1,numel(rho_values)+1);
    labels{1} = 'Full-state loop';
    for rho_index = 1:numel(rho_values)
        labels{rho_index+1} = ['$\rho=' ...
            num2str(rho_values(rho_index),'%.0e') '$'];
    end
end

function export_task5_figures(figures,output_dir)
    export_dir = fullfile(output_dir,'plots_eps');
    if ~exist(export_dir,'dir')
        mkdir(export_dir);
    end

    drawnow;
    previous_warning = warning('query', ...
        'MATLAB:print:ContentTypeImageSuggested');
    warning('off','MATLAB:print:ContentTypeImageSuggested');
    cleanup_warning = onCleanup(@() warning(previous_warning.state, ...
        'MATLAB:print:ContentTypeImageSuggested'));

    for figure_index = 1:numel(figures)
        fig = figures(figure_index);
        base_name = regexprep(fig.Name,'[^a-zA-Z0-9]+','_');
        base_name = regexprep(base_name,'^_+|_+$','');
        eps_path = fullfile(export_dir,[base_name '.eps']);
        png_path = fullfile(export_dir,[base_name '.png']);
        try
            exportgraphics(fig,eps_path, ...
                'ContentType','vector','BackgroundColor','w');
            exportgraphics(fig,png_path, ...
                'Resolution',300,'BackgroundColor','w');
        catch
            print(fig,fullfile(export_dir,base_name), ...
                '-depsc','-r300');
            print(fig,png_path,'-dpng','-r300');
        end
    end
end

function figures = Plot_Task6_Results(results)
% Build the comprehensive five-scenario comparison.

    task6 = results.task6;
    n_controllers = numel(task6.controller_labels);
    colors = turbo(n_controllers);
    line_styles = {'-','--','-.',':','-','--','-.',':'};
    time_label = '$t~[\mathrm{s}]$';
    figures = gobjects(1,3);

    figures(1) = figure('Color','w', ...
        'Name','Task 6 Comprehensive Metrics', ...
        'Position',[40 25 1450 900]);
    metric_values = {task6.rms_tracking_error, ...
        task6.maximum_control_effort,task6.final_parameter_error};
    metric_labels = {'$e_{\mathrm{RMS}}$', ...
        '$\max_t\|u(t)\|_2$', ...
        '$\|\widetilde{\Theta}(T)\|_F$'};
    metric_titles = {'RMS tracking error','Maximum control effort', ...
        'Final uncertainty-estimation error'};
    for metric_index = 1:3
        subplot(3,1,metric_index);
        values = max(metric_values{metric_index},1e-12);
        bars = bar(values,'grouped');
        for controller_index = 1:n_controllers
            bars(controller_index).FaceColor = colors(controller_index,:);
        end
        set(gca,'YScale','log','XTick',1:numel(task6.scenario_labels), ...
            'XTickLabel',task6.scenario_labels, ...
            'TickLabelInterpreter','latex');
        latex_style(metric_labels{metric_index},'', ...
            metric_titles{metric_index});
        if metric_index == 1
            legend(task6.controller_labels,'Interpreter','latex', ...
                'Location','northoutside','NumColumns',4,'FontSize',8);
        end
    end
    latex_super_title('Task 6: performance over all required scenarios');

    figures(2) = figure('Color','w', ...
        'Name','Task 6 Scenario Tracking Comparison', ...
        'Position',[60 30 1400 850]);
    for scenario_index = 1:numel(task6.scenarios)
        subplot(3,2,scenario_index); hold on;
        cases = task6.scenarios(scenario_index).cases;
        for controller_index = 1:n_controllers
            item = cases{controller_index};
            semilogy(item.t,max(vecnorm(item.e,2,2),1e-10), ...
                'Color',colors(controller_index,:), ...
                'LineStyle',line_styles{controller_index}, ...
                'LineWidth',1.25);
        end
        event_time = task6.scenarios(scenario_index).event_time;
        if isfinite(event_time)
            xline(event_time,'k--','LineWidth',1.0);
        end
        latex_style('$\|e(t)\|_2$',time_label, ...
            task6.scenario_labels{scenario_index});
        if scenario_index == 1
            legend(task6.controller_labels,'Interpreter','latex', ...
                'Location','best','NumColumns',2,'FontSize',7);
        end
    end
    latex_super_title('Task 6: tracking-error histories');

    step_index = find(strcmp(task6.scenario_labels,'Parameter step'),1);
    step_cases = task6.scenarios(step_index).cases;
    step_time = task6.scenarios(step_index).parameter_step_time;
    figures(3) = figure('Color','w', ...
        'Name','Task 6 Parameter Step Adaptation', ...
        'Position',[80 120 1400 520]);
    subplot(1,2,1); hold on;
    for controller_index = 1:n_controllers
        item = step_cases{controller_index};
        semilogy(item.t,max(vecnorm(item.e,2,2),1e-10), ...
            'Color',colors(controller_index,:), ...
            'LineStyle',line_styles{controller_index}, ...
            'LineWidth',1.35);
    end
    xline(step_time,'k--','Parameter change', ...
        'Interpreter','latex','LineWidth',1.0);
    latex_style('$\|e(t)\|_2$',time_label,'Tracking recovery');
    legend(task6.controller_labels,'Interpreter','latex', ...
        'Location','best','NumColumns',2,'FontSize',8);

    subplot(1,2,2); hold on;
    for controller_index = 1:n_controllers
        item = step_cases{controller_index};
        semilogy(item.t,max(item.Theta_error_norm,1e-10), ...
            'Color',colors(controller_index,:), ...
            'LineStyle',line_styles{controller_index}, ...
            'LineWidth',1.35);
    end
    xline(step_time,'k--','Parameter change', ...
        'Interpreter','latex','LineWidth',1.0);
    latex_style('$\|\widetilde{\Theta}(t)\|_F$', ...
        time_label,'Parameter re-identification');
    latex_super_title('Task 6: response to a sudden plant-parameter change');

    if isfield(results,'output_dir') && ~isempty(results.output_dir)
        export_task6_figures(figures,results.output_dir);
    end
end

function export_task6_figures(figures,output_dir)
    export_dir = fullfile(output_dir,'plots_eps');
    if ~exist(export_dir,'dir')
        mkdir(export_dir);
    end

    drawnow;
    previous_warning = warning('query', ...
        'MATLAB:print:ContentTypeImageSuggested');
    warning('off','MATLAB:print:ContentTypeImageSuggested');
    cleanup_warning = onCleanup(@() warning(previous_warning.state, ...
        'MATLAB:print:ContentTypeImageSuggested'));

    for figure_index = 1:numel(figures)
        fig = figures(figure_index);
        base_name = regexprep(fig.Name,'[^a-zA-Z0-9]+','_');
        base_name = regexprep(base_name,'^_+|_+$','');
        eps_path = fullfile(export_dir,[base_name '.eps']);
        png_path = fullfile(export_dir,[base_name '.png']);
        try
            exportgraphics(fig,eps_path, ...
                'ContentType','vector','BackgroundColor','w');
            exportgraphics(fig,png_path, ...
                'Resolution',300,'BackgroundColor','w');
        catch
            print(fig,fullfile(export_dir,base_name), ...
                '-depsc','-r300');
            print(fig,png_path,'-dpng','-r300');
        end
    end
end
