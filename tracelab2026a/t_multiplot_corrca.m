function figs = t_multiplot_corrca(cordat, cfg)
% Combines topoplots, boxplot, and timecourse into a composite figure per group
%--------------------------------------------------------------------------
% Configuration Parameters for t_multiplot_corrca(cfg)
%
% Optional (defaults listed):

%   cfg.components            → Which components to plot, e.g. [1 2 3] or 'sig_comps'
%                               If 'sig_comps', selects components with p-values < cfg.alpha
%                               (default: 1:3)
%   cfg.norm_method           → Method to normalize topoplot weights:
%                                 'zscore' → mean-centered & standardized
%                                 'minmax' → scaled to [0,1] range
%                               (default: 'zscore')
%   cfg.intensity_weighting   → If true, scales topoplot weights by ISC strength (default: true)
%   cfg.dot_style             → Landmark display style on topoplot:
%                                 'light'   → faint grey dots
%                                 'numbers' → numbered landmarks
%                                 'none'    → hide landmarks
%                                 <color>   → any valid color code
%                               (default: 'light')
%   cfg.visible               → Toggle figure display during generation (default: true)
%   cfg.colormap              → Colormap for topoplots (e.g. jet, parula) (default: jet)
%   cfg.c_label               → Colorbar label text for topoplots (default: 'Weight (a.u.)')
%   cfg.ylim                  → Color limits for topoplots:
%                                 [] = auto-range
%                                 [min max] = fixed range
%                                 'global' = use cfg.caxis
%                               (default: [])
%   cfg.symmetric             → If true, uses symmetric color range around zero (default: true)
%   cfg.output_format         → Format(s) for saving output figures:
%                                 e.g. 'png', 'pdf', 'svg', 'fig', etc.
%                               Can be string or cell array (default: {'png'})
%   cfg.output_dir            → Directory to save exported figures
%                               If empty, no saving performed (default: '')
%   cfg.alpha                 → Significance threshold (used if components = 'sig_comps') (default: 0.05)

% Notes:
% - Generates composite plots per group: topoplots, ISC boxplots, and timecourses.
% - Uses face landmark template from face_template.mat to render topographic maps.
% - Requires precomputed CorrCA results with fields: AvgA, ISC, ISC_persubject, ISC_persecond, etc.
%--------------------------------------------------------------------------

% --- Config defaults ---
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'components'), cfg.components = 1:3; end
if ~isfield(cfg, 'norm_method'), cfg.norm_method = 'zscore'; end
if ~isfield(cfg, 'dot_style'), cfg.dot_style = 'light'; end
if ~isfield(cfg, 'intensity_weighting'), cfg.intensity_weighting = true; end
if ~isfield(cfg, 'visible'), cfg.visible = true; end
if ~isfield(cfg, 'colormap'), cfg.colormap = jet; end
if ~isfield(cfg, 'c_label') || isempty(cfg.c_label)
    cfg.c_label = generate_colorbar_label(cfg);
end
if ~isfield(cfg, 'ylim'), cfg.ylim = []; end
if ~isfield(cfg, 'symmetric'), cfg.symmetric = true; end
if ~isfield(cfg, 'output_format'), cfg.output_format = {'png'}; end
if ~isfield(cfg, 'output_dir'), cfg.output_dir = ''; end
if ischar(cfg.output_format), cfg.output_format = {cfg.output_format}; end

vis_str = 'off'; if cfg.visible, vis_str = 'on'; end

% --- Load face template landmarks ---
fpath = fullfile(fileparts(mfilename('fullpath')), 'face_template.mat');
S = load(fpath); meanFace = S.meanFace;

groupNames = fieldnames(cordat);
figs = struct();

for i = 1:numel(groupNames)
    group = groupNames{i};
    res = cordat.(group);
    comps = cfg.components;
    comps = comps(comps <= size(res.AvgA,2));
    nComp = numel(comps);
    nCols = nComp + 1;

    % --- Normalize and scale weights ---
    weights_all = cell(1, nComp);
    all_weights = [];
    for j = 1:nComp
        w = res.AvgA(:, comps(j));
        switch lower(cfg.norm_method)
            case 'zscore', w = (w - mean(w)) / std(w);
            case 'minmax', w = (w - min(w)) / (max(w) - min(w) + eps);
        end
        if cfg.intensity_weighting && isfield(res, 'ISC')
            w = w * (res.ISC(comps(j)) / max(res.ISC(:)));
        end
        weights_all{j} = w;
        all_weights = [all_weights; w(:)];
    end

    % --- Determine clim ---
    if ischar(cfg.ylim) && strcmpi(cfg.ylim, 'global') && isfield(cfg, 'caxis')
        clim = cfg.caxis;
    elseif isnumeric(cfg.ylim) && numel(cfg.ylim) == 2
        clim = cfg.ylim;
    elseif isfield(cfg, 'symmetric') && cfg.symmetric
        lim = max(abs(all_weights(:)));
        clim = [-lim lim];
    else
        clim = [min(all_weights), max(all_weights)];
    end

    % --- Layout setup ---
    fig = figure('Color','w', 'Visible', vis_str, ...
        'Units','normalized', 'Position',[0.1 0.1 0.9 0.85]);
    t = tiledlayout(2, nCols, 'TileSpacing','compact', 'Padding','compact');

    % --- Topoplots: row 1, cols 1:nComp ---
    ax_topo = gobjects(1, nComp);
    for j = 1:nComp
        ax = nexttile(t, j);
        axes(ax);
        resTemp = res;
        resTemp.AvgA = weights_all{j};
        plot_corrca_topoplot(weights_all{j}, resTemp, comps(j), meanFace, cfg);
        caxis(clim);
        title(ax, sprintf('Comp %d', comps(j)), 'FontSize', 9, 'FontWeight','normal');
        ax_topo(j) = ax;
    end

    % --- Colorbar on middle topoplot ---
    cb = colorbar(ax_topo(ceil(nComp/2)), 'Location', 'southoutside');
    cb.Label.String = cfg.c_label;
    cb.Label.FontSize = 9;
    cb.Limits = clim;

    % --- Boxplot: row 1, col nCols ---
    axBox = nexttile(t, nCols);
    plot_corrca_boxplot(res, comps, cfg);
    title(axBox, 'ISC Boxplot');

    % --- Timecourse: row 2, spanning all columns ---
    axTime = nexttile(t, [1 nCols]);
    plot_corrca_timecourse(res, comps, cfg);
    title(axTime, sprintf('ISC Timecourse — Components %s', num2str(comps)));

    % --- Title and export ---
    sgtitle(fig, strrep(group, '_','\_'), 'FontWeight','bold');
    figs.(group) = fig;

    if ~isempty(cfg.output_dir)
        if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
        for f = 1:numel(cfg.output_format)
            fmt = lower(cfg.output_format{f});
            fname = sprintf('%s_combined.%s', group, fmt);
            exportgraphics(fig, fullfile(cfg.output_dir, fname), 'Resolution', 300);
        end
    end
end
end
