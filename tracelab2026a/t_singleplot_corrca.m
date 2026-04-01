function plot_archive = t_singleplot_corrca(cordata, cfg)

% This toolbox includes notBoxPlot.m by Rob Campbell (BSD license)

% =========================================================================
% cfg (Configuration Options)
% -------------------------------------------------------------------------
% cfg.components           : Components to plot (e.g. [1 2 3] or 'sig_comps')
% cfg.show                 : Cell array of plot types to display:
%                            {'topoplot','boxplot','timecourse'}
% cfg.norm_method          : Normalization method for weights:
%                            'none' (default), 'zscore', 'minmax'
% cfg.ylim                 : Color axis limits for topoplots:
%                            - Numeric [min max] for fixed scaling
%                            - 'global' to compute across all groups
% cfg.symmetric            : Whether to force symmetric color limits (true/false)
% cfg.dot_style            : Landmark overlay style:
%                            - 'light' = gray dots
%                            - 'numbers' = label each landmark
%                            - 'none' = no overlay
%                            - or any color (e.g. 'k' for black dots)
% cfg.intensity_weighting  : Apply ISC-based scaling to component weights (true/false)
% cfg.showTitle            : Add titles to figures (true/false)
% cfg.output_dir           : Folder path to save exported figures (if non-empty)
% cfg.output_format        : Cell array of formats to export:
%                            e.g. {'png'}, {'tiff','pdf'}, etc.
% cfg.colormap             : Colormap used for topoplots (e.g. 'jet', 'parula', custom Nx3 matrix)
% cfg.visible              : show plots in Matlab?
% =========================================================================

% --- Set defaults ---
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'components'),             cfg.components = 1:3; end
if ~isfield(cfg, 'show'),                   cfg.show = {'topoplot','boxplot','timecourse'}; end
if ~isfield(cfg, 'norm_method'),            cfg.norm_method = 'none'; end
if ~isfield(cfg, 'ylim'),                   cfg.ylim = []; end
if ~isfield(cfg, 'symmetric'),              cfg.symmetric = true; end
if ~isfield(cfg, 'dot_style'),              cfg.dot_style = 'light'; end
if ~isfield(cfg, 'intensity_weighting'),    cfg.intensity_weighting = false; end
if ~isfield(cfg, 'showTitle'),              cfg.showTitle = true; end
if ~isfield(cfg, 'output_dir'),             cfg.output_dir = ''; end
if ~isfield(cfg, 'output_format'),          cfg.output_format = {'png'}; end
if ischar(cfg.output_format),               cfg.output_format = {cfg.output_format}; end
if ~isfield(cfg, 'visible'), cfg.visible = true; end
if ~isfield(cfg, 'c_label') || isempty(cfg.c_label)
    cfg.c_label = generate_colorbar_label(cfg);
end

% --- Load face template landmarks ---
fpath = fullfile(fileparts(mfilename('fullpath')), 'face_template.mat');
S = load(fpath);
meanFace = S.meanFace;

groupNames = fieldnames(cordata);
plot_archive = struct();


% --- Global boxplot Y-limits across all groups ---
if ~isfield(cfg, 'boxplot_ylim') || isempty(cfg.boxplot_ylim)
    all_vals = [];
    for g = 1:numel(groupNames)
        res = cordata.(groupNames{g});
        comps = resolve_comps(res, cfg);
        comps = comps(comps <= size(res.ISC_persubject,1));
        all_vals = [all_vals; res.ISC_persubject(comps, :)'];
    end
    pad = 0.05 * range(all_vals(:));
    cfg.boxplot_ylim = [min(all_vals(:))-pad, max(all_vals(:))+pad];
end

% --- Global timecourse Y-limits across all groups ---
if ~isfield(cfg, 'timecourse_ylim') || isempty(cfg.timecourse_ylim)
    all_time_vals = [];
    for g = 1:numel(groupNames)
        res = cordata.(groupNames{g});
        comps = resolve_comps(res, cfg);
        comps = comps(comps <= size(res.ISC_persecond, 1));
        all_time_vals = [all_time_vals; res.ISC_persecond(comps, :)'];
    end
    pad = 0.05 * range(all_time_vals(:));
    cfg.timecourse_ylim = [min(all_time_vals(:)) - pad, max(all_time_vals(:)) + pad];
end



% --- Per-group loop ---
for g = 1:numel(groupNames)
    group = groupNames{g};
    res = cordata.(group);
    comps = resolve_comps(res, cfg);
    comps = comps(comps <= size(res.AvgA,2));
    entry = struct('topoplots', [], 'boxplot', [], 'timecourse', []);

    % --- Topoplot ---
    if ismember('topoplot', cfg.show)
        for c = 1:numel(comps)
            compIdx = comps(c);
            weights = res.AvgA(:, compIdx);

            switch lower(cfg.norm_method)
                case 'zscore'
                    weights = (weights - mean(weights)) / std(weights);
                case 'minmax'
                    weights = (weights - min(weights)) / (max(weights) - min(weights) + eps);
            end

            if cfg.intensity_weighting && isfield(res, 'ISC') && numel(res.ISC) >= compIdx
                scaleFactor = res.ISC(compIdx) / max(res.ISC(:));
                weights = weights * scaleFactor;
            end

            % --- Determine clim ---
            if ischar(cfg.ylim) && strcmpi(cfg.ylim, 'global') && isfield(cfg, 'caxis')
                clim = cfg.caxis;
            elseif isnumeric(cfg.ylim) && numel(cfg.ylim) == 2
                clim = cfg.ylim;
            else
                if cfg.symmetric
                    lim = max(abs(weights(:)));
                    clim = [-lim, lim];
                else
                    clim = [min(weights(:)), max(weights(:))];
                end
            end

            fig = figure('Color','w', 'Visible', cfg.visible);
            plot_corrca_topoplot(weights, res, compIdx, meanFace, cfg);

            cb = colorbar;
            cb.Label.String = cfg.c_label;
            cb.Label.FontSize = 9;
            cb.Limits = clim;

            if cfg.showTitle
                title(sprintf('%s — Component %d', strrep(group,'_',' '), compIdx), 'FontWeight','bold');
            end

            entry.topoplots{c} = fig;

            if ~isempty(cfg.output_dir)
                if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
                for f = 1:numel(cfg.output_format)
                    fmt = lower(cfg.output_format{f});
                    fname = sprintf('%s_comp%d_topoplot.%s', group, compIdx, fmt);
                    fpath = fullfile(cfg.output_dir, fname);
                    export_figure(fig, fpath, fmt);
                end
            end
        end
    end

    % --- Boxplot ---
    if ismember('boxplot', cfg.show)
        fig_w = max(400, 150 * numel(comps));  % Dynamic width based on num components
fig = figure('Color','w', 'Visible', cfg.visible, 'Position', [100 100 fig_w 300]);
        plot_corrca_boxplot(res, comps, cfg);
        if cfg.showTitle
            title(sprintf('%s — ISC Boxplot (Comps %s)', strrep(group,'_',' '), num2str(comps)), 'FontWeight','bold');
        end
        entry.boxplot = fig;
        if ~isempty(cfg.output_dir)
            for f = 1:numel(cfg.output_format)
                fmt = lower(cfg.output_format{f});
                fname = sprintf('%s_boxplot.%s', group, fmt);
                fpath = fullfile(cfg.output_dir, fname);
                export_figure(fig, fpath, fmt);
            end
        end
    end

    % --- Timecourse ---
    if ismember('timecourse', cfg.show)
        fig = figure('Color','w', 'Visible',cfg.visible);
        plot_corrca_timecourse(res, comps, cfg);
        if cfg.showTitle
            title(sprintf('%s — ISC Over Time (Comp %d)', strrep(group,'_',' '), comps(1)), 'FontWeight','bold');
        end
        entry.timecourse = fig;
        if ~isempty(cfg.output_dir)
            for f = 1:numel(cfg.output_format)
                fmt = lower(cfg.output_format{f});
                fname = sprintf('%s_timecourse.%s', group, fmt);
                fpath = fullfile(cfg.output_dir, fname);
                export_figure(fig, fpath, fmt);
            end
        end
    end

    plot_archive.(group) = entry;
end
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
% generate_colorbar_label - Create descriptive colorbar label based on settings

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