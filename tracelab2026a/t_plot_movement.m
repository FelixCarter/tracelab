function t_plot_movement(summary_data, cfg)
% T_PLOT_MOVEMENT - Visualize head movement data across participants and conditions
% Takes output from t_movement_summary as input
%
% cfg struct fields:
% ------------------
%   cfg.avg_groups      - Boolean. If true, average across groups/conditions on single plot. 
%                         If false, plot conditions separately. Default = false
%   cfg.tiles_per_plot  - Integer. Number of conditions to show per plot when plotting individually. 
%                         Default = 8
%   cfg.truncate_names  - Boolean. If true, truncate condition names to 15 characters. 
%                         Default = true
%   cfg.colormap        - String or Nx3 matrix. Colormap name ('lines', 'hsv', 'hot', etc.) 
%                         or custom RGB colors. Default = 'lines'
%   cfg.ylim            - [min max] vector, 'global', or 'off'. Y-axis limits. 
%                         'global' uses mean ± SEM range across all conditions. 
%                         'off' uses automatic per-plot limits. Default = 'global'
%   cfg.show_individual - Boolean. If true, show individual participant traces in light gray.
%                         Only applies when cfg.avg_groups = false. Default = false
%   cfg.error_display   - String. 'shaded' for filled area or 'whisker' for error bars.
%                         Default = 'shaded'
%
% Output:
% -------
%   Creates plots based on the configuration. No return value.
%
% Usage Examples:
% ---------------
%   % Basic individual condition plots
%   cfg_plot.avg_groups = false;
%   cfg_plot.tiles_per_plot = 6;
%   t_plot_movement(summary_data, cfg_plot);
%
%   % All conditions on one plot with custom colors
%   cfg_plot.avg_groups = true;
%   cfg_plot.colormap = 'hsv';
%   t_plot_movement(summary_data, cfg_plot);
%
%   % Bar chart comparison with custom y-limits
%   cfg_summary.avg_ppt = true;
%   summary_data = t_movement_summary(dat, cfg_summary);
%   cfg_plot.avg_groups = true;
%   cfg_plot.ylim = [0 5];
%   t_plot_movement(summary_data, cfg_plot);

    % Set default configuration
    if nargin < 2
        cfg = struct();
    end
    
    % Default values
    if ~isfield(cfg, 'avg_groups'), cfg.avg_groups = false; end
    if ~isfield(cfg, 'tiles_per_plot'), cfg.tiles_per_plot = 8; end
    if ~isfield(cfg, 'truncate_names'), cfg.truncate_names = true; end
    if ~isfield(cfg, 'colormap'), cfg.colormap = 'lines'; end
    if ~isfield(cfg, 'ylim'), cfg.ylim = 'global'; end
    if ~isfield(cfg, 'show_individual'), cfg.show_individual = false; end
    if ~isfield(cfg, 'error_display'), cfg.error_display = 'shaded'; end
    
    % Extract configuration from summary data
    data_cfg = summary_data.cfg;
    condition_names = summary_data.condition_names;
    
    % Determine plot type based on configurations
    if logical(data_cfg.avg_ppt) && logical(cfg.avg_groups)
        plot_bar_chart_conditions(summary_data, cfg);
    elseif ~logical(data_cfg.avg_ppt) && logical(cfg.avg_groups)
        plot_timecourse_all_conditions(summary_data, cfg);
    else
        plot_timecourse_individual(summary_data, cfg);
    end
end

function plot_timecourse_individual(summary_data, cfg)
% Plot individual timecourse for each condition
    
    condition_names = summary_data.condition_names;
    num_conditions = length(condition_names);
    
    % Calculate global y-limits based on MEAN data only (not individual participants)
    if ischar(cfg.ylim) && strcmp(cfg.ylim, 'global')
        all_mean_data = [];
        all_sem_data = [];
        
        for cond_idx = 1:num_conditions
            cond_name = condition_names{cond_idx};
            condition_data = summary_data.(cond_name);
            
            if isfield(condition_data, 'averaged_table') && ~isempty(condition_data.averaged_table)
                avg_table = condition_data.averaged_table;
                % Only use mean ± SEM for global range calculation
                all_mean_data = [all_mean_data; avg_table.MeanMovement];
                all_sem_data = [all_sem_data; avg_table.SEM];
            else
                % Calculate mean from participant data if no averaged table
                participant_tables = condition_data.participant_tables;
                has_valid_data = condition_data.has_valid_data;
                if any(has_valid_data)
                    valid_tables = participant_tables(has_valid_data);
                    
                    % Find common time points and calculate mean
                    all_times = [];
                    for ppt_idx = 1:length(valid_tables)
                        if ~isempty(valid_tables{ppt_idx})
                            all_times = union(all_times, valid_tables{ppt_idx}.Time);
                        end
                    end
                    
                    if ~isempty(all_times)
                        mean_movement = zeros(size(all_times));
                        for time_idx = 1:length(all_times)
                            values = [];
                            for ppt_idx = 1:length(valid_tables)
                                if ~isempty(valid_tables{ppt_idx})
                                    match_idx = abs(valid_tables{ppt_idx}.Time - all_times(time_idx)) < 1e-10;
                                    if any(match_idx)
                                        values = [values; valid_tables{ppt_idx}.Movement(match_idx)];
                                    end
                                end
                            end
                            mean_movement(time_idx) = mean(values, 'omitnan');
                        end
                        all_mean_data = [all_mean_data; mean_movement];
                    end
                end
            end
        end
        
        if ~isempty(all_mean_data)
            y_min = min(all_mean_data, [], 'omitnan');
            y_max = max(all_mean_data, [], 'omitnan');
            y_range = y_max - y_min;
            % Add 20% padding to see error bars properly
            global_ylim = [max(0, y_min - 0.2*y_range), y_max + 0.2*y_range];
        else
            global_ylim = [0 1]; % Fallback
        end
    elseif isnumeric(cfg.ylim) && length(cfg.ylim) == 2
        global_ylim = cfg.ylim;
    else
        global_ylim = []; % Auto per plot
    end
    
    % Calculate number of plots needed
    tiles_per_plot = cfg.tiles_per_plot;
    num_plots = ceil(num_conditions / tiles_per_plot);
    
    for plot_idx = 1:num_plots
        figure;
        
        % Calculate conditions for this plot
        start_cond = (plot_idx - 1) * tiles_per_plot + 1;
        end_cond = min(plot_idx * tiles_per_plot, num_conditions);
        conditions_this_plot = start_cond:end_cond;
        num_conditions_this_plot = length(conditions_this_plot);
        
        % Use 0.5N x 2 layout as requested
        rows = max(1, ceil(num_conditions_this_plot / 2));
        cols = min(2, num_conditions_this_plot);
        
        for local_idx = 1:num_conditions_this_plot
            cond_idx = conditions_this_plot(local_idx);
            cond_name = condition_names{cond_idx};
            condition_data = summary_data.(cond_name);
            
            subplot(rows, cols, local_idx);
            hold on;
            
            % Check if we have averaged data
            if isfield(condition_data, 'averaged_table') && ~isempty(condition_data.averaged_table)
                avg_table = condition_data.averaged_table;
                
                % Plot individual participant traces if requested
                if cfg.show_individual
                    participant_tables = condition_data.participant_tables;
                    has_valid_data = condition_data.has_valid_data;
                    if any(has_valid_data)
                        valid_tables = participant_tables(has_valid_data);
                        for ppt_idx = 1:length(valid_tables)
                            if ~isempty(valid_tables{ppt_idx})
                                ppt_table = valid_tables{ppt_idx};
                                plot(ppt_table.Time, ppt_table.Movement, 'Color', [0.85 0.85 0.85], ...
                                    'LineWidth', 0.5, 'HandleVisibility', 'off');
                            end
                        end
                    end
                end
                
                % Plot mean with error display
                if strcmp(cfg.error_display, 'shaded')
                    % Shaded error region
                    upper_bound = avg_table.MeanMovement + avg_table.SEM;
                    lower_bound = avg_table.MeanMovement - avg_table.SEM;
                    fill([avg_table.Time; flipud(avg_table.Time)], ...
                         [upper_bound; flipud(lower_bound)], ...
                         'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    % Mean line
                    plot(avg_table.Time, avg_table.MeanMovement, 'b-', 'LineWidth', 2);
                else
                    % Error bars
                    errorbar(avg_table.Time, avg_table.MeanMovement, avg_table.SEM, ...
                        'b-', 'LineWidth', 2, 'CapSize', 5);
                end
                
            else
                % Fallback: calculate mean from participant data
                participant_tables = condition_data.participant_tables;
                has_valid_data = condition_data.has_valid_data;
                
                if any(has_valid_data)
                    valid_tables = participant_tables(has_valid_data);
                    
                    % Plot individual participant traces if requested
                    if cfg.show_individual
                        for ppt_idx = 1:length(valid_tables)
                            if ~isempty(valid_tables{ppt_idx})
                                ppt_table = valid_tables{ppt_idx};
                                plot(ppt_table.Time, ppt_table.Movement, 'Color', [0.85 0.85 0.85], ...
                                    'LineWidth', 0.5, 'HandleVisibility', 'off');
                            end
                        end
                    end
                    
                    % Calculate and plot mean
                    all_times = [];
                    for ppt_idx = 1:length(valid_tables)
                        if ~isempty(valid_tables{ppt_idx})
                            all_times = union(all_times, valid_tables{ppt_idx}.Time);
                        end
                    end
                    
                    if ~isempty(all_times)
                        mean_movement = zeros(size(all_times));
                        sem_movement = zeros(size(all_times));
                        
                        for time_idx = 1:length(all_times)
                            values = [];
                            for ppt_idx = 1:length(valid_tables)
                                if ~isempty(valid_tables{ppt_idx})
                                    match_idx = abs(valid_tables{ppt_idx}.Time - all_times(time_idx)) < 1e-10;
                                    if any(match_idx)
                                        values = [values; valid_tables{ppt_idx}.Movement(match_idx)];
                                    end
                                end
                            end
                            mean_movement(time_idx) = mean(values, 'omitnan');
                            sem_movement(time_idx) = std(values, 'omitnan') / sqrt(length(values));
                        end
                        
                        % Plot with error display
                        if strcmp(cfg.error_display, 'shaded')
                            upper_bound = mean_movement + sem_movement;
                            lower_bound = mean_movement - sem_movement;
                            fill([all_times; flipud(all_times)], ...
                                 [upper_bound; flipud(lower_bound)], ...
                                 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                            plot(all_times, mean_movement, 'b-', 'LineWidth', 2);
                        else
                            errorbar(all_times, mean_movement, sem_movement, ...
                                'b-', 'LineWidth', 2, 'CapSize', 5);
                        end
                    end
                end
            end
            
            % Set y-axis limits
            if ~isempty(global_ylim)
                ylim(global_ylim);
            end

            % Set x-axis limits to match actual data range
            if isfield(condition_data, 'averaged_table') && ~isempty(condition_data.averaged_table)
                % Use averaged table if available
                xlim([min(condition_data.averaged_table.Time), max(condition_data.averaged_table.Time)]);
            elseif any(condition_data.has_valid_data)
                % Fallback: use first valid participant's time range
                valid_tables = condition_data.participant_tables(condition_data.has_valid_data);
                for t = 1:length(valid_tables)
                    if ~isempty(valid_tables{t})
                        first_valid = valid_tables{t};
                        xlim([min(first_valid.Time), max(first_valid.Time)]);
                        break;
                    end
                end
            end                        
            % Format title
            if cfg.truncate_names && length(cond_name) > 15
                display_name = [cond_name(1:12) '...'];
            else
                display_name = cond_name;
            end
            title(display_name, 'Interpreter', 'none');
            xlabel('Time (s)');
            ylabel('Nose energy');
            grid on;
            hold off;
        end
        
    end
end

function plot_timecourse_all_conditions(summary_data, cfg)
% Plot all conditions on single timecourse plot
    
    condition_names = summary_data.condition_names;
    num_conditions = length(condition_names);
    
    figure; hold on;
    
    % Define colors for different conditions
    if ischar(cfg.colormap)
        colors = feval(cfg.colormap, num_conditions);
    else
        colors = cfg.colormap;
    end
    
    legend_entries = {};
    valid_conditions = 0;
    
    for cond_idx = 1:num_conditions
        cond_name = condition_names{cond_idx};
        condition_data = summary_data.(cond_name);
        
        % Check if we have data to plot
        has_data = false;
        
        if isfield(condition_data, 'averaged_table') && ~isempty(condition_data.averaged_table)
            % Use averaged data
            avg_table = condition_data.averaged_table;
            x_data = avg_table.Time;
            y_data = avg_table.MeanMovement;
            err_data = avg_table.SEM;
            has_data = true;
            
        elseif any(condition_data.has_valid_data)
            % Calculate mean from participant data
            valid_tables = condition_data.participant_tables(condition_data.has_valid_data);
            
            if ~isempty(valid_tables)
                % Find common time points
                all_times = [];
                for ppt_idx = 1:length(valid_tables)
                    if ~isempty(valid_tables{ppt_idx})
                        all_times = union(all_times, valid_tables{ppt_idx}.Time);
                    end
                end
                
                if ~isempty(all_times)
                    mean_movement = zeros(size(all_times));
                    sem_movement = zeros(size(all_times));
                    
                    for time_idx = 1:length(all_times)
                        values = [];
                        for ppt_idx = 1:length(valid_tables)
                            if ~isempty(valid_tables{ppt_idx})
                                match_idx = abs(valid_tables{ppt_idx}.Time - all_times(time_idx)) < 1e-10;
                                if any(match_idx)
                                    values = [values; valid_tables{ppt_idx}.Movement(match_idx)];
                                end
                            end
                        end
                        mean_movement(time_idx) = mean(values, 'omitnan');
                        sem_movement(time_idx) = std(values, 'omitnan') / sqrt(length(values));
                    end
                    
                    x_data = all_times;
                    y_data = mean_movement;
                    err_data = sem_movement;
                    has_data = true;
                end
            end
        end
        
        if has_data
            valid_conditions = valid_conditions + 1;
            
            % Plot mean with error bars (SEM)
            errorbar(x_data, y_data, err_data, ...
                'Color', colors(cond_idx, :), 'LineWidth', 2, ...
                'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', colors(cond_idx, :));
            
            % Format legend entry
            if cfg.truncate_names && length(cond_name) > 15
                legend_entry = [cond_name(1:12) '...'];
            else
                legend_entry = cond_name;
            end
            legend_entries{end+1} = legend_entry;
        end
    end
    
    if valid_conditions > 0
        xlabel('Time (s)');
        ylabel('Nose energy');
        title('Head Movement Time Course - All Conditions');
        legend(legend_entries, 'Interpreter', 'none', 'Location', 'best');
        grid on;
    else
        title('No valid data to plot');
    end
    hold off;
end

function plot_bar_chart_conditions(summary_data, cfg)
% Plot bar chart showing average movement for each condition
% ERROR BARS NOW SHOW BETWEEN-PARTICIPANT VARIABILITY
    
    condition_names = summary_data.condition_names;
    num_conditions = length(condition_names);
    
    condition_means = zeros(num_conditions, 1);
    condition_sems = zeros(num_conditions, 1);
    valid_conditions = false(num_conditions, 1);
    
    for cond_idx = 1:num_conditions
        cond_name = condition_names{cond_idx};
        condition_data = summary_data.(cond_name);
        
        % Calculate participant-level overall means for between-participant variability
        participant_tables = condition_data.participant_tables;
        has_valid_data = condition_data.has_valid_data;
        
        if any(has_valid_data)
            valid_tables = participant_tables(has_valid_data);
            participant_overall_means = zeros(length(valid_tables), 1);
            
            for ppt_idx = 1:length(valid_tables)
                if ~isempty(valid_tables{ppt_idx})
                    % Calculate mean across time for each participant
                    participant_overall_means(ppt_idx) = mean(valid_tables{ppt_idx}.Movement, 'omitnan');
                else
                    participant_overall_means(ppt_idx) = NaN;
                end
            end
            
            % Remove NaN participants
            valid_participant_means = participant_overall_means(~isnan(participant_overall_means));
            
            if ~isempty(valid_participant_means) && length(valid_participant_means) > 1
                % Condition mean = mean of participant means (between-participant)
                condition_means(cond_idx) = mean(valid_participant_means, 'omitnan');
                % SEM = variability across participants (between-participant)
                condition_sems(cond_idx) = std(valid_participant_means, 'omitnan') / sqrt(length(valid_participant_means));
                valid_conditions(cond_idx) = true;
            else
                condition_means(cond_idx) = NaN;
                condition_sems(cond_idx) = NaN;
                valid_conditions(cond_idx) = false;
            end
        else
            condition_means(cond_idx) = NaN;
            condition_sems(cond_idx) = NaN;
            valid_conditions(cond_idx) = false;
        end
    end
    
    figure;
    
    if any(valid_conditions)
        % Plot only valid conditions
        valid_idx = find(valid_conditions);
        num_valid = length(valid_idx);
        
        % Get colors for bars
        if ischar(cfg.colormap)
            colors = feval(cfg.colormap, num_valid);
        else
            colors = cfg.colormap;
            % If custom colors provided but not enough for all conditions, cycle through them
            if size(colors, 1) < num_valid
                color_indices = mod(0:num_valid-1, size(colors, 1)) + 1;
                colors = colors(color_indices, :);
            end
        end
        
        % Create colored bars - NO LEGEND
        if num_valid == 1
            % Single bar case
            bar(valid_idx, condition_means(valid_conditions), 'FaceColor', colors(1,:));
        else
            % Multiple bars case - simplified without legend handles
            bar_handles = bar(valid_idx, condition_means(valid_conditions), 'FaceColor', 'flat');
            for i = 1:num_valid
                bar_handles.CData(i,:) = colors(i,:);
            end
        end
        hold on;
        
        % Add error bars (now showing between-participant variability)
        errorbar(valid_idx, condition_means(valid_conditions), condition_sems(valid_conditions), ...
            'k.', 'LineWidth', 1.5, 'CapSize', 10);
        
        % Format x-axis labels
        x_labels = condition_names(valid_conditions);
        if cfg.truncate_names
            for i = 1:length(x_labels)
                if length(x_labels{i}) > 15
                    x_labels{i} = [x_labels{i}(1:12) '...'];
                end
            end
        end
        
        xticks(valid_idx);
        xticklabels(x_labels);
        xtickangle(45);
        
        ylabel('Nose energy');
        title('Average Head Movement by Condition');
        grid on;
        
        % Set y-axis limits if specified
        if isnumeric(cfg.ylim) && length(cfg.ylim) == 2
            ylim(cfg.ylim);
        end
    else
        text(0.5, 0.5, 'No averaged data available', 'HorizontalAlignment', 'center', ...
            'Units', 'normalized', 'FontSize', 12);
        title('No Data to Plot');
    end
    hold off;
end