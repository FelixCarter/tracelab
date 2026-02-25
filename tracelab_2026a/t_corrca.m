function results = t_corrca(dat, cfg)
% t_corrca - Run CorrCA-based ISC per group with optional significance testing
%--------------------------------------------------------------------------
% Supported landmark transformations:
%   cfg.combine_landmarks → 'pca_z', 'amplitude', 'vecnorm' (default: 'vecnorm')
% Configuration options (defaults listed):
%   cfg.norm_landmarks     → Mean-center landmark coordinates (default: true)
%   cfg.gamma              → Regularization parameter (default: 0.1)
%   cfg.n_sec              → Window length in seconds (default: 3)
%   cfg.fs                 → Sampling rate in Hz (default: inferred)
%   cfg.pca_scale          → Scale mode for PCA-based compression (default: 'none')
%   cfg.zscore             → Z-score mode: 'none', 'pre', 'post', or 'both' (default: 'none')
%   cfg.rescale            → Whether to min-max normalize features to [0,1] 
%                            prior to covariance computation (default: true)
%                            Ensures all landmarks contribute equally
%--------------------------------------------------------------------------

if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'norm_landmarks'), cfg.norm_landmarks = true; end
if ~isfield(cfg, 'gamma'), cfg.gamma = 0.1; end
if ~isfield(cfg, 'n_sec'), cfg.n_sec = 3; end
if ~isfield(cfg, 'combine_landmarks'), cfg.combine_landmarks = 'vecnorm'; end
if ~isfield(cfg, 'pca_scale'), cfg.pca_scale = 'none'; end
if ~isfield(cfg, 'zscore'), cfg.zscore = 'none'; end
if ~isfield(cfg, 'rescale'), cfg.rescale = true; end

valid_methods = {'pca_z', 'amplitude', 'vecnorm', 'kpm'};
if ~ismember(cfg.combine_landmarks, valid_methods)
    error('Unsupported preprocessing method: %s. Choose one of: %s', ...
        cfg.combine_landmarks, strjoin(valid_methods, ', '));
end

valid_zscore_modes = {'none', 'pre', 'post', 'both'};
if ~ismember(cfg.zscore, valid_zscore_modes)
    error('Invalid cfg.zscore mode: "%s". Choose from: %s', ...
        cfg.zscore, strjoin(valid_zscore_modes, ', '));
end

groups = fieldnames(dat);
results = struct();

for g = 1:numel(groups)
    group = groups{g};
    entries = dat.(group);
    N = numel(entries);

    % Verify time alignment across participants
    times = cellfun(@(x) x.XY.time, entries, 'UniformOutput', false);
    aligned = all(cellfun(@(t) isequal(t, times{1}), times));
    if ~aligned
        warning('Skipping group "%s": inconsistent time vectors.', group);
        continue;
    end

    T = length(times{1});
    D = width(entries{1}.XY) - 1;
    X = nan(T, D, N);

    for i = 1:N
        XY = table2array(entries{i}.XY(:, 2:end));
        if cfg.norm_landmarks
            XY = XY - mean(XY, 1);
        end
        X(:, :, i) = double(XY);
    end

    % ⚖️ Pre-transform z-score normalization
    if ismember(cfg.zscore, {'pre', 'both'})
        for d = 1:D
            for n = 1:N
                X(:, d, n) = zscore(X(:, d, n));
            end
        end
    end

    % Transform landmarks
    norm_cfg = struct('type', cfg.pca_scale);
    X = transform_landmarks(X, cfg.combine_landmarks, norm_cfg);
    D = size(X, 2);

    % ⚖️ Post-transform z-score normalization
    if ismember(cfg.zscore, {'post', 'both'})
        for l = 1:D
            for n = 1:N
                X(:, l, n) = zscore(X(:, l, n));
            end
        end
    end

    % Infer sampling rate from metadata
    if ~isfield(cfg, 'fs')
        cfg.fs = entries{1}.Settings.target_fs;
        fprintf('Using fs=%.2f Hz from Settings for "%s"\n', cfg.fs, group);
    end

    try
        [ISC, ISC_persubject, ISC_persecond, W, A] = isc_tracelab(X, cfg);

        % Initialize results structure
        results.(group) = struct( ...
            'ISC', ISC, ...
            'ISC_persubject', ISC_persubject, ...
            'ISC_persecond', ISC_persecond, ...
            'W', W, ...
            'A', A, ...
            'AvgW', W, ...
            'AvgA', A, ...
            'ParticipantIDs', {cellfun(@(x) x.ParticipantID, entries, 'UniformOutput', false)} ...
        );

        % Save runtime configuration info for traceability
        settings = struct( ...
            'fs', cfg.fs, ...
            'n_sec', cfg.n_sec, ...
            'combine_landmarks', cfg.combine_landmarks, ...
            'norm_landmarks', cfg.norm_landmarks, ...
            'zscore_mode', cfg.zscore, ...
            'gamma', cfg.gamma ...  % ADD THIS LINE
        );
        if isfield(cfg, 'pca_scale')
            settings.pca_scale = cfg.pca_scale;
        end
        results.(group).Settings = settings;

        if isfield(cfg, 'pca_scale')
            settings.pca_scale = cfg.pca_scale;
        end
        results.(group).Settings = settings;

        fprintf('CorrCA complete for "%s" (%d participants, %d features)\n', group, N, D);

    catch ME
        fprintf('Failed on group "%s": %s\n', group, ME.message);
    end
end
end

function Xtrans = transform_landmarks(Xraw, method, norm_cfg)
% TRANSFORM_LANDMARKS - Compress or reduce x-y landmark pairs using specified method
%
% Inputs:
%   Xraw     → T × D × N array (time × features × subjects), where D = 136
%   method   → string method: 'none' | 'pca_z' | 'amplitude' | 'vecnorm'
%   norm_cfg → (optional) struct with field:
%                .type = 'zscore' | 'unit' | 'none'  (default: 'none')
%
% Output:
%   Xtrans   → T × D′ × N transformed matrix (e.g. D′ = 68 for compressed modes)

[T, D, N] = size(Xraw);
num_landmarks = D / 2;
Xtrans = [];
if nargin < 3 || ~isfield(norm_cfg, 'type'), norm_cfg.type = 'none'; end

switch lower(method)
    case 'none'
        Xtrans = Xraw;

    case 'pca_z'
        Xtrans = nan(T, num_landmarks, N);
        for n = 1:N
            for l = 1:num_landmarks
                x = Xraw(:, 2*l-1, n);
                y = Xraw(:, 2*l, n);

                % Optional normalization
                switch lower(norm_cfg.type)
                    case 'zscore'
                        x = (x - mean(x)) / std(x); y = (y - mean(y)) / std(y);
                    case 'unit'
                        mag = sqrt(x.^2 + y.^2);
                        x = x ./ max(mag, eps); y = y ./ max(mag, eps);
                end

                pair = [x, y];
                if any(isnan(pair(:)))
                    proj = nan(T,1);
                else
                    [~, proj, ~] = pca(pair, 'NumComponents', 1);
                end
                Xtrans(:, l, n) = proj;
            end
        end

    case 'vecnorm'
        Xtrans = nan(T, num_landmarks, N);
        for n = 1:N
            for l = 1:num_landmarks
                x = Xraw(:, 2*l-1, n);
                y = Xraw(:, 2*l, n);
                Xtrans(:, l, n) = sqrt(x.^2 + y.^2);
            end
        end

        case 'amplitude'
            Xtrans = nan(T - 1, num_landmarks, N);  % 💥 Shorten time dimension
            for n = 1:N
                for l = 1:num_landmarks
                    x = Xraw(:, 2*l-1, n);
                    y = Xraw(:, 2*l, n);
                    dx = diff(x);  % no padding
                    dy = diff(y);
                    Xtrans(:, l, n) = sqrt(dx.^2 + dy.^2);
                end
            end

            case 'kpm'
            Xtrans = nan(T, num_landmarks, N);
            for n = 1:N
                neutral = mean(Xraw(:, :, n), 1);  % [1 × D]
                for l = 1:num_landmarks
                    x = Xraw(:, 2*l-1, n) - neutral(2*l-1);
                    y = Xraw(:, 2*l, n)   - neutral(2*l);
                    Xtrans(:, l, n) = sqrt(x.^2 + y.^2);  % Euclidean displacement
                end
            end


    otherwise
        error('Unknown transformation method: %s', method);
end
end

function [ISC, ISC_persubject, ISC_persecond, W, A] = isc_tracelab(X, cfg)
% isc_tracelab - CorrCA-based inter-subject correlation analysis
%
% Inputs:
%   X   : [T × D × N] data matrix (e.g. mean-centered landmark or EEG data)
%   cfg : struct with optional fields:
%         - gamma: regularization (default 0.1)
%         - n_sec: time window in seconds (default 3)
%         - fs   : sampling frequency (default 5)
%
% Outputs:
%   ISC               : eigenvalues (ISC per component)
%   ISC_persubject    : ISC per component × subject
%   ISC_persecond     : ISC per second per component
%   W, A              : spatial filters and forward models

% ------------------ Config ------------------
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'gamma'), cfg.gamma = 0.1; end
if ~isfield(cfg, 'n_sec'), cfg.n_sec = 3; end
if ~isfield(cfg, 'fs'),    cfg.fs = 5; end

gamma = cfg.gamma;
Nsec  = cfg.n_sec;
fs    = cfg.fs;

[T, D, N] = size(X);
if cfg.rescale
    X = rescale(X);  
end

% ------------------ Cross-covariances Rij ------------------
X2 = reshape(X, T, D*N);  % ensure 2D flattening
C = cov(X2);

expectedSize = D * N;
if size(C,1) ~= expectedSize
    error('Unexpected cov size: expected [%d × %d], got [%d × %d]', ...
        expectedSize, expectedSize, size(C));
end

Rij = permute(reshape(C, [D N D N]), [1 3 2 4]);  % D × D × N × N

Rw = 1/N * sum(Rij(:, :, 1:N+1:N*N), 3);
Rb = 1/(N*(N-1)) * (sum(Rij(:,:,:), 3) - N * Rw);
Rw_reg = (1 - gamma) * Rw + gamma * mean(eig(Rw)) * eye(D);
%Rw_reg = (1 - gamma) * Rw + gamma * eye(D);

[W, ISC] = eig(Rb, Rw_reg);
[ISC, idx] = sort(diag(ISC), 'descend');
W = W(:, idx);
A = Rw * W / (W' * Rw * W);

% ------------------ ISC per subject ------------------
ISC_persubject = zeros(length(ISC), N);
for i = 1:N
    Rw_i = 0; Rb_i = 0;
    for j = 1:N
        if i == j, continue; end
        Rw_i = Rw_i + 1/(N-1) * (Rij(:,:,i,i) + Rij(:,:,j,j));
        Rb_i = Rb_i + 1/(N-1) * (Rij(:,:,i,j) + Rij(:,:,j,i));
    end
    ISC_persubject(:,i) = diag(W' * Rb_i * W) ./ diag(W' * Rw_i * W);
end

% ------------------ ISC per second ------------------
ISC_persecond = [];
winSize = Nsec * fs;
fprintf("cfg.n_sec = %.2f\n", cfg.n_sec);
fprintf("Smoothing window length = %.2f\n", winSize);
for t = 1:floor((T - winSize)/fs)
    idx = (1:winSize) + (t-1)*fs;
    Xt = X(idx, :, :);
    Xt2 = reshape(Xt, winSize, D*N);
    Ct = cov(Xt2);
    if size(Ct,1) ~= expectedSize, continue; end
    Rij_t = permute(reshape(Ct, [D N D N]), [1 3 2 4]);
    Rw_t = 1/N * sum(Rij_t(:, :, 1:N+1:N*N), 3);
    Rb_t = 1/(N*(N-1)) * (sum(Rij_t(:,:,:), 3) - N * Rw_t);
    ISC_persecond(:,t) = diag(W' * Rb_t * W) ./ diag(W' * Rw_t * W);
end
end

