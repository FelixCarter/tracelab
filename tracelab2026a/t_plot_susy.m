function fig_handles = t_plot_susy(summary_out, cfg)
% t_plot_susy - Visualize SUSY synchrony over time OR as bar plots for aggregated data

    % t_plot_susy - Visualize SUSY synchrony over time with error bars or shaded regions
    % takes output from t_susy_summarise as input
    
    % cfg struct fields:
    % ------------------
    %   cfg.measure          - String. Measure to plot, e.g. 'Z' or 'Z_abs'. Default = 'Z'
    %   cfg.compare_pseudo   - Boolean. If true, plot matching pseudo measure
    %   (e.g. 'Z_pseudo_abs'). Default = true
    %   cfg.combine_segments - Integer ≥ 1. Number of segments to combine when binning. Default = 1
    %   cfg.error_type       - String. 'sem', 'std', or 'ci'. Determines error calculation. Default = 'sem'
    %   cfg.error_display    - String. 'whisker' or 'shaded'. How to display variability. Default = 'shaded'
    %   cfg.avg_level        - String. 'dyad' or 'participant'. Level to average at. Default = 'dyad'
    %   cfg.segment          - Required. Segment length in seconds
    %   cfg.fs               - Required. Sampling rate in Hz
    %   cfg.output_dir       - String. Directory to save plots. Default = '' (i.e. don't save)
    %   cfg.output_format    - String or cell array. File format(s) to save as: 'png', 'pdf', 'fig', 'tiff', 'jpg', etc. Default = 'png'
    %   cfg.ylim             - Either [min max] or 'global' for consistent y-limits across plots. Default = auto
    %   cfg.x_format         - String. 'seconds' or 'minutes'. Controls x-axis label format. Default = 'seconds'
    %   cfg.legend_location  - String. Where to put legend? Default = 'northeast'
    %   cfg.title            - Boolean. If set to false, plots won't have titles

if ~isfield(cfg, 'measure'),           cfg.measure = 'Z'; end
if ~isfield(cfg, 'compare_pseudo'),    cfg.compare_pseudo = true; end
if ~isfield(cfg, 'combine_segments'),  cfg.combine_segments = 1; end
if ~isfield(cfg, 'error_type'),        cfg.error_type = 'sem'; end
if ~isfield(cfg, 'error_display'),     cfg.error_display = 'shaded'; end
if ~isfield(cfg, 'avg_level'),         cfg.avg_level = 'dyad'; end
if ~isfield(cfg, 'output_dir'),        cfg.output_dir = ''; end
if ~isfield(cfg, 'output_format'),     cfg.output_format = 'png'; end
if ~isfield(cfg, 'ylim'),              cfg.ylim = []; end
if ~isfield(cfg, 'x_format'),          cfg.x_format = 'seconds'; end
if ~isfield(cfg, 'legend_location'),   cfg.legend_location = 'northeast'; end
if ~isfield(cfg, 'title'),             cfg.title = true; end
if ~isfield(cfg, 'plot_type'),         cfg.plot_type = 'auto'; end 
if ~isfield(cfg, 'segment'),           cfg.segment = 15; end  

% ===== DETERMINE PLOT TYPE =====
groups = fieldnames(summary_out);
is_timecourse = false;
has_dyad_columns = false;

for g = 1:numel(groups)
    T = summary_out.(groups{g});
    if ~isempty(T)
        % Check for timecourse data
        if ismember('Segment', T.Properties.VariableNames)
            is_timecourse = true;
        end
        
        % Check for dyad columns
        if all(ismember({'Dyad1', 'Dyad2'}, T.Properties.VariableNames))
            has_dyad_columns = true;
        end
        break;
    end
end

% Auto-detect plot type if not specified
if strcmp(cfg.plot_type, 'auto')
    if is_timecourse && has_dyad_columns
        cfg.plot_type = 'timecourse';
    else
        cfg.plot_type = 'bar';
    end
end

% ===== BRANCH TO APPROPRIATE PLOTTING FUNCTION =====
if strcmp(cfg.plot_type, 'bar')
    fig_handles = plot_aggregated_bars(summary_out, cfg);
else
    fig_handles = plot_timecourse(summary_out, cfg, is_timecourse, has_dyad_columns);
end
end

%% ===== TIME COURSE PLOTTING =====
function fig_handles = plot_timecourse(summary_out, cfg, is_timecourse, has_dyad_columns)

groups = fieldnames(summary_out);
fig_handles = cell(numel(groups), 1);

% Global y-limit computation (modified)
globalY = [];
if ischar(cfg.ylim) && strcmpi(cfg.ylim, 'global')
    all_vals = []; all_errs = [];
    for g = 1:numel(groups)
        T = summary_out.(groups{g});
        if isempty(T), continue; end
        
        % Handle different data formats
        if is_timecourse
            if ismember('Segment', T.Properties.VariableNames)
                T.Bin = floor((T.Segment - 1) / cfg.combine_segments) + 1;
            else
                T.Bin = ones(height(T), 1);
            end
        else
            % For aggregated data, just use single point
            T.Bin = ones(height(T), 1);
        end
        
        if has_dyad_columns && strcmpi(cfg.avg_level, 'participant')
            T1 = T; T1.Unit = T1.Dyad1;
            T2 = T; T2.Unit = T2.Dyad2;
            T = [T1; T2];
        elseif has_dyad_columns
            T.Unit = strcat(T.Dyad1, "-", T.Dyad2);
        else
            % If no dyad columns, use Participant or create dummy units
            if ismember('Participant', T.Properties.VariableNames)
                T.Unit = T.Participant;
            else
                T.Unit = arrayfun(@(x) sprintf('Unit%d', x), 1:height(T), 'UniformOutput', false)';
            end
        end
        
        bins = unique(T.Bin);
        [m, e] = compute_trace(T, bins, cfg.measure, cfg, has_dyad_columns);
        all_vals = [all_vals; m(:)];
        all_errs = [all_errs; e(:)];
    end
    if all(isfinite(all_vals)) && all(isfinite(all_errs))
        globalY = [min(all_vals - all_errs), max(all_vals + all_errs)];
    end
end

% Plotting loop
for g = 1:numel(groups)
    group = groups{g};
    T = summary_out.(group);
    if isempty(T), warning("Group '%s' is empty", group); continue; end
    if ~ismember(cfg.measure, T.Properties.VariableNames)
        error("Measure '%s' not found.", cfg.measure);
    end
    
    % Prepare data structure based on available columns
    if has_dyad_columns && strcmpi(cfg.avg_level, 'participant')
        T1 = T; T1.Unit = T1.Dyad1;
        T2 = T; T2.Unit = T2.Dyad2;
        T = [T1; T2];
    elseif has_dyad_columns
        T.Unit = strcat(T.Dyad1, "-", T.Dyad2);
    elseif ismember('Participant', T.Properties.VariableNames)
        T.Unit = T.Participant;
    else
        T.Unit = arrayfun(@(x) sprintf('Unit%d', x), 1:height(T), 'UniformOutput', false)';
    end
    
    % Handle binning
    if is_timecourse && ismember('Segment', T.Properties.VariableNames)
        T.Bin = floor((T.Segment - 1) / cfg.combine_segments) + 1;
        bins = unique(T.Bin);
        dur = cfg.segment;  
        t = ((bins - 1) + 0.5) * cfg.combine_segments * dur;
    else
        % Single point for aggregated data
        T.Bin = ones(height(T), 1);
        bins = 1;
        t = 0.5;  % Center at 0.5
    end
    
    [mean_vals, err_vals] = compute_trace(T, bins, cfg.measure, cfg, has_dyad_columns);
    
    % Create figure
    if is_timecourse
        fig_name = sprintf('Synchrony – %s', group);
    else
        fig_name = sprintf('Aggregated Synchrony – %s', group);
    end
    fig = figure('Name', fig_name, 'Color', 'w'); hold on;
    
    % Plot with appropriate styling
    if length(t) > 1
        % Timecourse plot
        h_real = plot_trace(t, mean_vals, err_vals, cfg.error_display, [0.2 0.4 0.8]);
        
        % Axis formatting for timecourse
        if is_timecourse
            tick_pos = (0:cfg.combine_segments:(max(bins))) * dur;
            set(gca, 'XTick', tick_pos);
            if numel(tick_pos) > 30
                skip = ceil(numel(tick_pos) / 30);
                label_idx = 1:skip:numel(tick_pos);
            else
                label_idx = 1:numel(tick_pos);
            end
            
            if strcmpi(cfg.x_format, 'minutes')
                all_labels = arrayfun(@(s) datestr(seconds(s), 'MM:SS'), tick_pos, 'UniformOutput', false);
            else
                all_labels = cellstr(string(tick_pos));
            end
            
            final_labels = repmat({''}, size(all_labels));
            final_labels(label_idx) = all_labels(label_idx);
            set(gca, 'XTickLabel', final_labels);
            
            if strcmpi(cfg.x_format, 'minutes')
                xlabel('Time (MM:SS)');
            else
                xlabel(sprintf('Time (s) [%.1fs × %d]', dur, cfg.combine_segments));
            end
            
            xlim([min(tick_pos), max(tick_pos)]);
        end
    else
        % Single point plot (scatter/bar hybrid)
        bar_width = 0.6;
        x_pos = 0.5;
        
        % Draw bar with error
        bar(x_pos, mean_vals, bar_width, 'FaceColor', [0.2 0.4 0.8], ...
            'EdgeColor', 'k', 'LineWidth', 1.5);
        errorbar(x_pos, mean_vals, err_vals, 'k.', 'LineWidth', 2, 'CapSize', 10);
        
        % Also add individual points
        if height(T) <= 50  % Only show individual points if not too many
            jitter = (rand(height(T), 1) - 0.5) * bar_width * 0.8;
            scatter(x_pos + jitter, T.(cfg.measure), 40, ...
                'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.8);
        end
        
        xlim([0, 1]);
        set(gca, 'XTick', x_pos, 'XTickLabel', {group});
        xlabel('Group');
    end
    
    % Y-axis label
    if strcmpi(cfg.measure, 'Z')
        label_core = 'Z: non-abs';
    else
        label_core = cfg.measure;
    end
    ylabel(sprintf('Synchrony (%s)', label_core));
    
    % Title
    if cfg.title
        if is_timecourse
            title_str = sprintf('Group: %s  |  Level: %s  |  Measure: %s', ...
                group, cfg.avg_level, label_core);
        else
            title_str = sprintf('Aggregated Synchrony: %s  |  Measure: %s', ...
                group, label_core);
        end
        title(strrep(title_str, '_', ' '), 'Interpreter', 'none');
    end
    
    % Zero line
    yline(0, '--k', 'LineWidth', 1);
    
    % Set y-limits
    if isnumeric(cfg.ylim) && ~isempty(cfg.ylim)
        ylim(cfg.ylim);
    elseif ~isempty(globalY)
        ylim(globalY);
    end
    
    box on; grid on;
    
    % Legend
    if length(t) > 1  % Only for timecourse
        err_label = struct('sem','± Std Error','std','± Std Dev','ci','± 95% CI');
        legend_entries = {['Real ' err_label.(lower(cfg.error_type))]};
        legend_handles = h_real;
        
        % Optional pseudo comparison
        if cfg.compare_pseudo
            pseudo_var = strrep(cfg.measure, 'Z', 'Z_pseudo');
            if ismember(pseudo_var, T.Properties.VariableNames)
                [p_vals, p_err] = compute_trace(T, bins, pseudo_var, cfg, has_dyad_columns);
                h_pseudo = plot_trace(t, p_vals, p_err, cfg.error_display, [0.8 0.2 0.2]);
                legend_entries{2} = ['Pseudo ' err_label.(lower(cfg.error_type))];
                legend_handles(2) = h_pseudo;
            end
        end
        
        legend(legend_handles, legend_entries, 'Location', cfg.legend_location);
    end
    
    fig_handles{g} = fig;
    
    % Save figure
    save_figure(fig, cfg, group, 'timecourse');
end

end

%% ===== AGGREGATED BAR PLOTS =====
function fig_handles = plot_aggregated_bars(summary_out, cfg)
% Plot aggregated data as bar charts with error bars

groups = fieldnames(summary_out);
fig_handles = cell(1, 1);  % Single figure for all groups

% Extract data
group_data = struct();
for g = 1:numel(groups)
    T = summary_out.(groups{g});
    if isempty(T), continue; end
    
    % Get real synchrony
    if ismember(cfg.measure, T.Properties.VariableNames)
        group_data(g).name = groups{g};
        group_data(g).real_vals = T.(cfg.measure);
        group_data(g).real_mean = mean(T.(cfg.measure), 'omitnan');
        group_data(g).real_err = get_err(T.(cfg.measure), cfg.error_type);
    end
    
    % Get pseudo synchrony if requested
    if cfg.compare_pseudo
        pseudo_var = strrep(cfg.measure, 'Z', 'Z_pseudo');
        if ismember(pseudo_var, T.Properties.VariableNames)
            group_data(g).pseudo_vals = T.(pseudo_var);
            group_data(g).pseudo_mean = mean(T.(pseudo_var), 'omitnan');
            group_data(g).pseudo_err = get_err(T.(pseudo_var), cfg.error_type);
        end
    end
    
    % Store participant IDs for labeling
    if ismember('Participant', T.Properties.VariableNames)
        group_data(g).participants = T.Participant;
    end
end

% Remove empty groups
group_data = group_data(arrayfun(@(x) ~isempty(x.name), group_data));

% Create figure
fig = figure('Name', 'SUSY Aggregated Results', 'Color', 'w', ...
    'Position', [100, 100, 800, 600]);
fig_handles{1} = fig;
hold on;

n_groups = length(group_data);
bar_width = 0.35;
x_pos = 1:n_groups;

% Plot bars
real_means = [group_data.real_mean];
real_errs = [group_data.real_err];

h_real = bar(x_pos - bar_width/2, real_means, bar_width, ...
    'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);

% Add error bars
errorbar(x_pos - bar_width/2, real_means, real_errs, 'k.', ...
    'LineWidth', 1.5, 'CapSize', 10);

% Plot pseudo synchrony if available
if cfg.compare_pseudo && all(isfield(group_data, {'pseudo_mean', 'pseudo_err'}))
    pseudo_means = [group_data.pseudo_mean];
    pseudo_errs = [group_data.pseudo_err];
    
    h_pseudo = bar(x_pos + bar_width/2, pseudo_means, bar_width, ...
        'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    errorbar(x_pos + bar_width/2, pseudo_means, pseudo_errs, 'k.', ...
        'LineWidth', 1.5, 'CapSize', 10);
end

% Add individual data points
for g = 1:n_groups
    % Jitter x-positions
    n_points = length(group_data(g).real_vals);
    jitter = (rand(n_points, 1) - 0.5) * bar_width * 0.8;
    
    % Plot real synchrony points
    scatter(x_pos(g) - bar_width/2 + jitter, group_data(g).real_vals, 40, ...
        'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.8);
    
    % Plot pseudo synchrony points if available
    if cfg.compare_pseudo && isfield(group_data(g), 'pseudo_vals')
        scatter(x_pos(g) + bar_width/2 + jitter, group_data(g).pseudo_vals, 40, ...
            'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.8);
    end
end

% Formatting
xticks(x_pos);
xticklabels(strrep({group_data.name}, '_', ' '));
xlabel('Condition');

% Y-axis label
if strcmpi(cfg.measure, 'Z')
    ylabel('Synchrony (Z)');
else
    ylabel(sprintf('Synchrony (%s)', cfg.measure));
end

% Title
if cfg.title
    title_str = sprintf('Aggregated SUSY Results (n=%d groups)', n_groups);
    if cfg.compare_pseudo
        title_str = [title_str, ' - Real vs Pseudo'];
    end
    title(strrep(title_str, '_', ' '), 'Interpreter', 'none');
end

% Zero line
yline(0, '--k', 'LineWidth', 1);

% Grid and box
box on; grid on;

% Legend
if cfg.compare_pseudo && exist('h_pseudo', 'var')
    legend([h_real, h_pseudo], {'Real Synchrony', 'Pseudo Synchrony'}, ...
        'Location', cfg.legend_location);
else
    legend(h_real, 'Real Synchrony', 'Location', cfg.legend_location);
end

% Set y-limits if specified
if ~isempty(cfg.ylim) && isnumeric(cfg.ylim)
    ylim(cfg.ylim);
else
    % Auto-adjust y-limits to include error bars
    all_vals = [];
    for g = 1:n_groups
        all_vals = [all_vals; group_data(g).real_vals];
        if cfg.compare_pseudo && isfield(group_data(g), 'pseudo_vals')
            all_vals = [all_vals; group_data(g).pseudo_vals];
        end
    end
    y_range = [min(all_vals), max(all_vals)];
    y_padding = diff(y_range) * 0.1;
    ylim([y_range(1)-y_padding, y_range(2)+y_padding]);
end

% Save figure
save_figure(fig, cfg, 'aggregated_results', 'bar');
end

%% ===== HELPER FUNCTIONS =====
function [mean_vals, err_vals] = compute_trace(T, bins, varname, cfg, has_dyad_columns)
    mean_vals = zeros(size(bins));
    err_vals  = zeros(size(bins));
    
    for i = 1:numel(bins)
        subT = T(T.Bin == bins(i), :);
        
        if has_dyad_columns && ismember('Unit', subT.Properties.VariableNames)
            % Use dyad/participant level averaging
            [~, ~, idx] = unique(subT.Unit);
            dyad_vals = accumarray(idx, subT.(varname), [], @(x) mean(x, 'omitnan'));
            vals_to_use = dyad_vals;
        else
            % Use all values directly
            vals_to_use = subT.(varname);
        end
        
        mean_vals(i) = mean(vals_to_use, 'omitnan');
        err_vals(i)  = get_err(vals_to_use, cfg.error_type);
    end
end

function e = get_err(vals, type)
    vals = vals(~isnan(vals));
    n = length(vals);
    if n < 2
        warning('Insufficient data points (n=%d) to compute error. Returning NaN.', n);
        e = NaN;
        return;
    end
    
    switch lower(type)
        case 'sem'
            e = std(vals) / sqrt(n);
        case 'std'
            e = std(vals);
        case 'ci'
            e = 1.96 * std(vals) / sqrt(n);
        otherwise
            e = std(vals);
    end
end

function h = plot_trace(t, mean_vals, err_vals, mode, color)
    switch lower(mode)
        case 'whisker'
            h = errorbar(t, mean_vals, err_vals, '-o', ...
                'Color', color, 'LineWidth', 1.8, ...
                'MarkerFaceColor', color);
        case 'shaded'
            hold on;
            fill_x = [t; flipud(t)];
            fill_y = [mean_vals + err_vals; flipud(mean_vals - err_vals)];
            fill(fill_x, fill_y, color, ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none');
            h = plot(t, mean_vals, '-o', ...
                'Color', color, 'LineWidth', 2, ...
                'MarkerFaceColor', color);
        otherwise
            error("Unknown error_display mode: %s", mode);
    end
end

function save_figure(fig, cfg, group_name, plot_type)
    if isempty(cfg.output_dir), return; end
    
    if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
    
    formats = cellstr(cfg.output_format);
    for f = 1:numel(formats)
        fmt = lower(formats{f});
        fname = sprintf('%s_%s_%s_%s.%s', ...
            group_name, plot_type, cfg.avg_level, cfg.measure, fmt);
        filepath = fullfile(cfg.output_dir, fname);
        
        switch fmt
            case 'png',  saveas(fig, filepath);
            case 'pdf',  print(fig, filepath, '-dpdf', '-bestfit');
            case {'tif','tiff'}, print(fig, filepath, '-dtiff');
            case 'fig',  savefig(fig, filepath);
            case {'jpg','jpeg'}, print(fig, filepath, '-djpeg');
            case 'eps',  print(fig, filepath, '-depsc');
            case 'svg',  print(fig, filepath, '-dsvg');
            otherwise,   warning("Unsupported format: %s", fmt);
        end
    end
end