function results = t_corrca(dat, cfg)
% t_corrca - Run CorrCA-based ISC per group with optional significance testing
%--------------------------------------------------------------------------
% Supported landmark transformations:
%   cfg.combine_landmarks → 'pca_z', 'amplitude', 'vecnorm' (default: 'vecnorm')
% Configuration options (defaults listed):
%   cfg.norm_landmarks        → Mean-center landmark coordinates (default: true)
%   cfg.gamma                 → Regularization parameter (default: 0.1)
%   cfg.n_sec                 → Window length in seconds (default: 3)
%   cfg.fs                    → Sampling rate in Hz (default: inferred)
%   cfg.pca_scale             → Scale mode for PCA-based compression (default: 'none')
%   cfg.zscore                → Z-score mode: 'none', 'pre', 'post', or 'both' (default: 'none')
%   cfg.sigtest               → Run permutation testing (default: false)
%   cfg.n_permutations        → Number of permutations if sigtest=true (default: 1000)
%   cfg.sig_components        → Components to test (default: 1:3)
%   cfg.alpha                 → Significance level (default: 0.05)
%   cfg.remove_head_movement  → Subtract nose tip after transformation to isolate local movements (default: true)
%--------------------------------------------------------------------------

if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'norm_landmarks'), cfg.norm_landmarks = true; end
if ~isfield(cfg, 'gamma'), cfg.gamma = 0.1; end
if ~isfield(cfg, 'n_sec'), cfg.n_sec = 3; end
if ~isfield(cfg, 'combine_landmarks'), cfg.combine_landmarks = 'vecnorm'; end
if ~isfield(cfg, 'pca_scale'), cfg.pca_scale = 'none'; end
if ~isfield(cfg, 'zscore'), cfg.zscore = 'none'; end
if ~isfield(cfg, 'remove_head_movement'), cfg.remove_head_movement = true; end

% Significance testing defaults
if ~isfield(cfg, 'sigtest'), cfg.sigtest = false; end
if ~isfield(cfg, 'n_permutations'), cfg.n_permutations = 1000; end
if ~isfield(cfg, 'sig_components'), cfg.sig_components = 1:3; end
if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end

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
    
    % Skip empty groups
    if N == 0
        warning('Skipping group "%s": no participants.', group);
        continue;
    end

    % Verify time alignment across participants
    times = cellfun(@(x) x.XY.time, entries, 'UniformOutput', false);
    aligned = all(cellfun(@(t) isequal(t, times{1}), times));
    if ~aligned
        warning('Skipping group "%s": inconsistent time vectors.', group);
        continue;
    end

    T = length(times{1});
    D = width(entries{1}.XY) - 1;  % Should be 136 (68 landmarks * 2 coordinates)
    X = nan(T, D, N);

    for i = 1:N
        XY = table2array(entries{i}.XY(:, 2:end));
        
        % NOTE: Head movement correction moved to AFTER landmark transformation
        % to avoid creating linear dependencies in raw XY coordinates
        
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

    % Transform landmarks (e.g., 136 channels -> 68 magnitude channels)
    norm_cfg = struct('type', cfg.pca_scale);
    X = transform_landmarks(X, cfg.combine_landmarks, norm_cfg);
    D = size(X, 2);  % Now D = 68 (number of landmarks)

    % ========== HEAD MOVEMENT CORRECTION (APPLIED AFTER TRANSFORMATION) ==========
    % Subtract nose tip movement from all landmarks to isolate local facial movements
    % Nose tip is landmark 34 (index 34 in the transformed 68-channel space)
    if cfg.remove_head_movement && D >= 34
        nose_ref = X(:, 34, :);  % T × 1 × N
        X = X - nose_ref;        % Subtract nose from all landmarks
        % This removes global head movement while preserving local expression differences
    end
    % ============================================================================

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

        % Store the preprocessed data for significance testing if needed
        if cfg.sigtest
            X_processed = rescale(X);  % Store the exact data used by isc_tracelab
        end

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
            'gamma', cfg.gamma, ...
            'remove_head_movement', cfg.remove_head_movement ...
        );
        if isfield(cfg, 'pca_scale')
            settings.pca_scale = cfg.pca_scale;
        end
        results.(group).Settings = settings;

        % Run significance testing if requested
        if cfg.sigtest
            fprintf('Running permutation tests for "%s" (%d permutations)...\n', group, cfg.n_permutations);
            
            sig_components = cfg.sig_components;
            sig_components = sig_components(sig_components <= size(W, 2));
            
            sig_results = struct();
            
            for comp_idx = 1:length(sig_components)
                comp = sig_components(comp_idx);
                fprintf('  Component %d: ', comp);
                
                % Project original data to get true correlation
                orig_projected = zeros(T, N);
                for subj = 1:N
                    orig_projected(:, subj) = squeeze(X_processed(:, :, subj)) * W(:, comp);
                end
                orig_corr_mat = corr(orig_projected);
                orig_corr = mean(orig_corr_mat(triu(true(N), 1)));
                
                % Permutation test
                null_dist = zeros(cfg.n_permutations, 1);
                
                for perm = 1:cfg.n_permutations
                    % Circular shift the PROCESSED data
                    X_shifted = nan(size(X_processed));
                    for subj = 1:N
                        shift = randi(T);
                        X_shifted(:, :, subj) = circshift(squeeze(X_processed(:, :, subj)), shift);
                    end
                    
                    % Project shifted data using original W
                    projected_shifted = zeros(T, N);
                    for subj = 1:N
                        projected_shifted(:, subj) = squeeze(X_shifted(:, :, subj)) * W(:, comp);
                    end
                    
                    % Compute correlation for permuted data
                    corr_mat_shifted = corr(projected_shifted);
                    null_dist(perm) = mean(corr_mat_shifted(triu(true(N), 1)));
                end
                
                % Calculate significance
                p_value = (sum(null_dist >= orig_corr) + 1) / (cfg.n_permutations + 1);
                significant = p_value < cfg.alpha;
                
                % Store significance results
                sig_results(comp).component_number = comp;
                sig_results(comp).ISC = ISC(comp);
                sig_results(comp).actual_correlation = orig_corr;
                sig_results(comp).p_value = p_value;
                sig_results(comp).significant = significant;
                sig_results(comp).null_distribution = null_dist;
                
                fprintf('ISC=%.4f, r=%.4f, p=%.4f, sig=%d\n', ...
                    ISC(comp), orig_corr, p_value, significant);
            end
            
            results.(group).significance = sig_results;
        end

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
%   method   → string method: 'none' | 'pca_z' | 'amplitude' | 'vecnorm' | 'kpm'
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
        Xtrans = nan(T - 1, num_landmarks, N);  % Shorten time dimension
        for n = 1:N
            for l = 1:num_landmarks
                x = Xraw(:, 2*l-1, n);
                y = Xraw(:, 2*l, n);
                dx = diff(x);
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
                Xtrans(:, l, n) = sqrt(x.^2 + y.^2);
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
X = rescale(X);  % unit variance across full matrix

% ------------------ Cross-covariances Rij ------------------
X2 = reshape(X, T, D*N);
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