function plot_corrca_topoplot(weights, res, compIdx, meanFace, cfg)
% Plot a topoplot for a single component using corrca data + config

% --- Normalize weights ---
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

% --- Color limits ---
if isfield(cfg, 'ylim') && ischar(cfg.ylim) && strcmpi(cfg.ylim, 'global') && isfield(cfg, 'caxis')
    clim = cfg.caxis;
elseif isfield(cfg, 'ylim') && isnumeric(cfg.ylim) && numel(cfg.ylim) == 2
    clim = cfg.ylim;
else
    % Per-component auto-scaling
    if isfield(cfg, 'symmetric') && cfg.symmetric
        lim = max(abs(weights(:)));
        clim = [-lim, lim];
    else
        clim = [min(weights(:)), max(weights(:))];
    end
end

% --- Interpolate ---
F = scatteredInterpolant(meanFace(:,1), meanFace(:,2), weights, 'natural', 'none');
[xq, yq] = meshgrid( ...
    linspace(min(meanFace(:,1))-10, max(meanFace(:,1))+10, 400), ...
    linspace(min(meanFace(:,2))-10, max(meanFace(:,2))+10, 400));
vq = F(xq, yq);

% --- Plot ---
contourf(xq, yq, vq, 40, 'LineColor', 'none');
if isfield(cfg, 'colormap') && ~isempty(cfg.colormap)
    colormap(cfg.colormap);
else
    colormap(jet);
end
caxis(clim);
axis equal ij off;
set(gca, 'XLim', [min(xq(:)), max(xq(:))]);
set(gca, 'YLim', [min(yq(:)), max(yq(:))]);

% --- Landmark overlay ---
hold on;
switch lower(cfg.dot_style)
    case 'light'
        scatter(meanFace(:,1), meanFace(:,2), 10, [0.75 0.75 0.75], 'filled');
    case 'numbers'
        for i = 1:size(meanFace,1)
            text(meanFace(i,1), meanFace(i,2), sprintf('%d', i), ...
                 'FontSize', 8, 'HorizontalAlignment','center', 'Color','k');
        end
    case 'none'
        % Do nothing
    otherwise
        scatter(meanFace(:,1), meanFace(:,2), 10, cfg.dot_style, 'filled');
end
hold off;

% --- Title ---
if isfield(cfg, 'showTitle') && cfg.showTitle
    title(sprintf('Component %d', compIdx), 'FontWeight','bold');
end
end