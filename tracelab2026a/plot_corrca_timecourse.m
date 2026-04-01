function plot_corrca_timecourse(res, comp_idx_vec, cfg)
% plot_corrca_timecourse - Timecourse of ISC for selected components

Y = res.ISC_persecond(comp_idx_vec, :)';  % [T × selected components]
% 🌐 Pull timing info from Settings
if isfield(res, 'Settings')
    fs = res.Settings.fs;
    n_sec = res.Settings.n_sec;
else
    fs = 1;     % fallback
    n_sec = 1;
end

stepSize_sec = 1;  % assuming step = fs frames → 1 second shift
t = (0:size(Y,1)-1) * stepSize_sec;

% Optional: center timestamps around window midpoint
% t = ((0:size(Y,1)-1) + 0.5) * n_sec;


% --- Determine Y-limits ---
if isfield(cfg, 'timecourse_ylim') && ~isempty(cfg.timecourse_ylim)
    y_lims = cfg.timecourse_ylim;
else
    pad = 0.05 * range(Y(:));
    y_lims = [min(Y(:))-pad, max(Y(:))+pad];
end

% --- Plot ---
ax = gca;
cla(ax); hold(ax, 'on');
plot(t, Y, 'LineWidth', 1.2);
yline(0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.75);

legend(ax, "Comp " + string(comp_idx_vec), 'Location','northeast', 'Box','off');
title(ax, sprintf('ISC Timecourse — Components %s', num2str(comp_idx_vec)));
xlabel(ax, 'Time (s)'); ylabel(ax, 'ISC');
ylim(ax, y_lims);
xlim(ax, [min(t), max(t)]);
xticks(ax, min(t):10:max(t));
set(ax, 'XMinorTick','on');
ax.XAxis.MinorTickValues = min(t):1:max(t);

% --- Gridlines like overview plot ---
y_vals = y_lims(1):0.01:y_lims(2);
grid_lines = gobjects(numel(y_vals),1);
for i = 1:numel(y_vals)
    y = y_vals(i);
    grid_lines(i) = yline(ax, y, '-', 'Color', [0.85 0.85 0.85], ...
                          'LineWidth', 0.5, 'HandleVisibility', 'off');
end
uistack(grid_lines, 'bottom');  % ensures grid is behind traces
end

function comps = resolve_comps(res, cfg)
% resolve_comps - Determine which components to plot

if ischar(cfg.components) && strcmpi(cfg.components, 'sig_comps')
    if isfield(res, 'ISC_stats') && isfield(res.ISC_stats, 'pValues')
        if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end
        p = res.ISC_stats.pValues;
        comps = find(p < cfg.alpha);
        if isempty(comps)
            warning('No significant components found — defaulting to Comp 1 only');
            comps = 1;
        end
    else
        error('No p-value information found in res.ISC_stats');
    end
else
    comps = cfg.components;
end
end

function export_figure(fig, fpath, fmt)
try
    % Prevent automatic flipping during export
    set(fig, 'InvertHardcopy', 'off');
    drawnow;  % ⏱️ Ensure rendering is complete before saving

    % Handle each format explicitly
    switch lower(fmt)
        case 'png'
            print(fig, fpath, '-dpng', '-r300');
        case 'pdf'
            print(fig, fpath, '-dpdf', '-bestfit');
        case {'tiff', 'tif'}
            print(fig, fpath, '-dtiff', '-r300');
        case {'jpg', 'jpeg'}
            print(fig, fpath, '-djpeg', '-r300');
        case 'eps'
            print(fig, fpath, '-depsc');
        case 'svg'
            print(fig, fpath, '-dsvg');
        case 'fig'
            savefig(fig, fpath);
        otherwise
            warning('Unsupported format: %s', fmt);
    end
catch ME
    warning(' Failed to save figure: %s\n%s', fpath, ME.message);
end
end

function label = generate_colorbar_label(cfg)
% generate_colorbar_label - Create specific label based on transformation settings
%
% INPUTS:
%   cfg - Configuration structure with norm_method, intensity_weighting, symmetric
%
% OUTPUTS:
%   label - Descriptive colorbar label

    % Base cases
    if strcmpi(cfg.norm_method, 'none') && ~cfg.intensity_weighting
        if isfield(cfg, 'symmetric') && cfg.symmetric
            label = 'Component weight (symmetric scaling)';
        else
            label = 'Component weight';
        end
        return;
    end
    
    % Build descriptive parts
    parts = {};
    
    % Add normalization method
    if strcmpi(cfg.norm_method, 'zscore')
        parts{end+1} = 'z-scored';
    elseif strcmpi(cfg.norm_method, 'minmax')
        parts{end+1} = 'min-max normalized';
    % 'none' case already handled above
    end
    
    % Add intensity weighting
    if isfield(cfg, 'intensity_weighting') && cfg.intensity_weighting
        parts{end+1} = 'ISC-weighted';
    end
    
    % Add symmetric scaling
    if isfield(cfg, 'symmetric') && cfg.symmetric
        parts{end+1} = 'symmetric';
    end
    
    % Construct final label
    if isempty(parts)
        label = 'Component weight';
    elseif length(parts) == 1
        label = ['Component weight (', parts{1}, ')'];
    else
        label = ['Component weight (', strjoin(parts, ', '), ')'];
    end
end