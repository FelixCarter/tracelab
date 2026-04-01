function [results, susy_extra] = t_susy(dat, cfg)
    % t_susy - Unified SUSY analysis (dyadic and mvSUSY modes)
    % based on the Surrogate Synchrony R package by Tschacher 2023
    % Tschacher W (2023). SUSY: Surrogate Synchrony. R package version 0.1.1, https://wtschacher.github.io/SUSY/.
    % credit to Prof. Tschacher for developing the SUSY procedure
    % This function takes output from t_preproc as input
%--------------------------------------------------------------------------
% Configuration Parameters for t_susy(cfg)
%
% Optional (defaults listed below):

%   cfg.datatype              → Name of signal field to analyze (e.g. 'HeadMovement' or emotion type)
%                               Must match field name in participant structs or Emotions table
%                               (default: 'HeadMovement')
%   cfg.segment               → Segment length in seconds used for SUSY computation (default: 3)
%   cfg.fs                    → Sampling rate in Hz for time series data (default: 1)
%   cfg.maxlag                → Maximum time lag (in seconds) considered for synchrony alignment (default: 1)
%   cfg.permutation           → Whether to include surrogate permutations for significance testing (default: true)
%   cfg.restrict_surrogates  → If true, limits surrogates to group-compatible pairings (default: false)
%   cfg.total_surrogates     → Number of surrogate pairings to generate when permutation testing is enabled (default: 500)
%   cfg.susy_mode             → Analysis mode:
%                                'dyadic' → pairwise synchrony comparisons (default)
%                                'mv'     → multivariate synchrony estimation via mvSUSY
% Notes:
% - Supports both direct signal fields and emotion-detection fields stored in `Emotions`.
% - Automatically aligns signal length across participants in each group.
% - Outputs a results table per group containing synchrony metrics and (optionally) permutation stats.
%--------------------------------------------------------------------------

    % Set default configuration
    default_cfg = struct( ...
        'datatype', 'HeadMovement', ...
        'segment', 3, ...
        'fs', 1, ...
        'maxlag', 1, ...
        'permutation', true, ...
        'restrict_surrogates', false, ...
        'total_surrogates', 500, ...
        'susy_mode', 'dyadic' ...
    );

    % Merge user config with defaults
    if nargin < 2 || isempty(cfg)
        cfg = default_cfg;
    else
        cfg = mergeStruct(default_cfg, cfg);
    end

       fprintf('=== STARTING t_susy (%s mode) ===\n', cfg.susy_mode);


    % Initialize results
    groups = fieldnames(dat);
    results = struct();
    susy_extra = struct();

    for g = 1:numel(groups)
        group = groups{g};
        fprintf('\n=== PROCESSING GROUP: %s ===\n', group);

        try
            % 1. Extract data
            participants = dat.(group);
            num_participants = numel(participants);
            first_part = participants{1};

            if isfield(first_part, 'Emotions') && ismember(cfg.datatype, first_part.Emotions.Properties.VariableNames)
                data_source = 'Emotions';
            elseif isfield(first_part, cfg.datatype)
                data_source = 'direct';
            else
                error('Data type "%s" not found for group "%s".', cfg.datatype, group);
            end

            max_length = getMaxLength(participants, data_source, cfg.datatype);
            X = NaN(max_length, num_participants);
            ids = strings(num_participants, 1);

            for p = 1:num_participants
                if strcmp(data_source, 'Emotions')
                    X(:,p) = participants{p}.Emotions.(cfg.datatype);
                else
                    X(:,p) = participants{p}.(cfg.datatype);
                end
                ids(p) = getParticipantID(participants{p}, p);
            end

            % 2. Run analysis based on mode
            switch lower(cfg.susy_mode)
                case 'dyadic'
                    fprintf('-- Running dyadic SUSY for %s --\n', group);
                    [outTable, summary_stats, ES_table] = susy_core(X, cfg, ids); 

                    % Store summary stats separately
                    if ~exist('all_summary_stats', 'var')
                        all_summary_stats = struct();
                    end
                    all_summary_stats.(group) = summary_stats;

                case 'mv'
                    fprintf('-- Running mvSUSY for %s --\n', group);
                    error('mvSUSY mode not yet implemented');
                    % Uncomment when implemented:
                    % mv_out = mv_susy_core(X, cfg);
                    % outTable = mv_susy_summary_table(mv_out);
                    % summary_stats = table();
                    % ES_table = table();

                otherwise
                    error('Unknown susy_mode: %s. Use "dyadic" or "mv".', cfg.susy_mode);
            end

            % Main output (backward compatible)
            results.(group) = outTable;
            
            % Extra output with summary and ES
            susy_extra.(group).summary = summary_stats;
            susy_extra.(group).ES_table = ES_table;

        catch ME
            fprintf('!! ERROR in group "%s": %s !!\n', group, ME.message);
            results.(group) = table();
            susy_extra.(group).summary = table();
            susy_extra.(group).ES_table = table();
        end
    end

    fprintf('\n=== ANALYSIS COMPLETE (%s mode) ===\n', cfg.susy_mode);
    assignin('base', 'susy_cfg', cfg);
end

function cfg = mergeStruct(default, user)
    fields = fieldnames(default);
    for i = 1:length(fields)
        if ~isfield(user, fields{i})
            user.(fields{i}) = default.(fields{i});
        end
    end
    cfg = user;
end

function max_len = getMaxLength(participants, source, datatype)
    max_len = 0;
    for p = 1:numel(participants)
        if strcmp(source, 'Emotions')
            len = height(participants{p}.Emotions);
        else
            len = length(participants{p}.(datatype));
        end
        max_len = max(max_len, len);
    end
end

function id = getParticipantID(participant, default_id)
    if isfield(participant, 'ParticipantID')
        id = string(participant.ParticipantID);
    else
        id = sprintf('Participant%d', default_id);
    end
end

function [tbl, summary_stats, ES_table] = susy_core(X, cfg, ids)
    % SUSY core analysis - PRECISE MATCH to R implementation
    % Returns:
    %   tbl: Detailed per-segment, per-lag table
    %   summary_stats: Summary statistics with ES values
    %   ES_table: ES and ES_pseudo values for each pair
    
    % cfg PARAMETERS:
    % cfg.fs: Sampling frequency (Hz) - converts time to samples
    % cfg.segment: Segment duration (seconds) - each analysis window
    % cfg.maxlag: Max lag for cross-correlation (seconds) - ±range
    % cfg.permutation: If true: all possible pairs, false: consecutive pairs
    % cfg.restrict_surrogates: Surrogate method (see below)
    % cfg.total_surrogates: Total # surrogates (used if restrict_surrogates=true)
    
    % restrict_surrogates modes:
    % false (default): Each segment vs ALL other segments (N×(N-1) comparisons)
    % true: Each segment vs fixed # of surrogates (total_surrogates/N per segment)
    
    fs = cfg.fs;
    range = round(cfg.segment * fs);
    maxlagHz = round(cfg.maxlag * fs);
    
    % WARNING CHECKS
    if cfg.maxlag == 0
        warning('maxlag = 0: ES will be Inf/NaN (standard deviation = 0 with only 1 lag point). Use maxlag >= 1 for meaningful ES.');
    elseif cfg.maxlag == 1
        warning('maxlag = 1: ES_lead1 and ES_lead2 will be Inf/NaN (standard deviation = 0 for single-lag leads). Use maxlag >= 3 for lead analysis.');
    end
    
    lags = -maxlagHz:maxlagHz;
    lagtimes2 = length(lags);

    % Clean columns
    validCols = ~all(isnan(X));
    X = X(:, validCols);
    ids = ids(validCols);

    % Generate dyad pairs (matches R's pairs logic exactly)
    if cfg.permutation
        pairs = nchoosek(1:size(X,2), 2);
    else
        pairs = reshape(1:size(X,2), 2, [])';
        if mod(size(X,2), 2) ~= 0
            error('When permutation=false, X must have even number of columns (R requirement)');
        end
    end

    % Preallocate summary structure (matches R's as.data.frame output)
    summary_stats = table();
    
    % Preallocate ES table
    ES_table = table();
    
    % Main detailed table
    tbl = table();

    for p = 1:size(pairs,1)
        x = X(:, pairs(p,1));
        y = X(:, pairs(p,2));

        valid = ~isnan(x) & ~isnan(y);
        x = x(valid);
        y = y(valid);
        len = length(x);
        size_val = len;  % matches R's 'size' variable

        numSeg = floor(len / range);
        % R uses: numberEpochen = round(size/range-0.499999)
        numberEpochen = round(len/range - 0.499999);
        if numberEpochen < 2
            warning('Skipping pair %d-%d: too short for segmenting.', pairs(p,1), pairs(p,2));
            continue
        end

        segX = reshape(x(1:numberEpochen*range), range, numberEpochen)';
        segY = reshape(y(1:numberEpochen*range), range, numberEpochen)';

        % Initialize as in R code
        meanccorrReal = zeros(1, lagtimes2);
        meanccorrPseudo = zeros(1, lagtimes2);
        meanccorrRealZ = zeros(1, lagtimes2);
        meanccorrRealZNotAbs = zeros(1, lagtimes2);
        meanccorrPseudoZ = zeros(1, lagtimes2);
        meanccorrPseudoZNotAbs = zeros(1, lagtimes2);
        nReal = zeros(numberEpochen, 1);
        nPseudo = zeros(numberEpochen, 1);

        % Variables for detailed output (per segment, per lag)
        Z_real_signed   = nan(numberEpochen, lagtimes2);
        Z_real_abs      = nan(numberEpochen, lagtimes2);
        Z_pseudo_signed = nan(numberEpochen, lagtimes2);
        Z_pseudo_abs    = nan(numberEpochen, lagtimes2);
        
        % NEW: Accumulators for segment-specific pseudo values
        pseudo_signed_accum = zeros(numberEpochen, lagtimes2);
        pseudo_abs_accum = zeros(numberEpochen, lagtimes2);
        pseudo_counts = zeros(numberEpochen, 1);

        % Determine surrogate method as in R
        if cfg.restrict_surrogates
            anzahlPseudosProEpoche = floor(cfg.total_surrogates / numberEpochen);
            if anzahlPseudosProEpoche < 1
                warning('Number of pseudos per total > number of segments, setting to 1');
                anzahlPseudosProEpoche = 1;
            end
            if anzahlPseudosProEpoche >= numberEpochen
                warning('surrogates.total > number of segments, setting to numberEpochen-1');
                anzahlPseudosProEpoche = numberEpochen - 1;
            end
            if anzahlPseudosProEpoche == numberEpochen - 1
                % R sets restrict.surrogates = FALSE in this case
                cfg.restrict_surrogates = false;
            end
        else
            anzahlPseudosProEpoche = numberEpochen - 1;
        end

        % Main loop over segments (Epochen)
        if cfg.restrict_surrogates
            % METHOD 1: restrict.surrogates = TRUE
            isOcc = zeros(numberEpochen, 1);
            n = 0;
            
            for s = 1:numberEpochen
                a = segX(s,:) - mean(segX(s,:));
                b = segY(s,:) - mean(segY(s,:));
                
                if std(a) < 1e-6 || std(b) < 1e-6
                    continue;
                end
                
                % Real correlation for this segment
                r_real = xcorr(a, b, maxlagHz, 'coeff');
                z_real = fisherZ(r_real);
                
                Z_real_signed(s,:) = z_real;
                Z_real_abs(s,:) = abs(z_real);
                
                % Update meanccorrRealZ (accumulate Z values as in R)
                for j = 1:lagtimes2
                    meanccorrRealZ(j) = meanccorrRealZ(j) + abs(z_real(j));
                    meanccorrRealZNotAbs(j) = meanccorrRealZNotAbs(j) + z_real(j);
                    meanccorrReal(j) = meanccorrReal(j) + abs(r_real(j));
                    nReal(s) = nReal(s) + abs(r_real(j));
                end
                nReal(s) = nReal(s) / lagtimes2;
                
                % Generate surrogates for this segment
                for h = 1:anzahlPseudosProEpoche
                    % Find unique surrogate segment (matches R's repeat loop)
                    while true
                        v = randi(numberEpochen);
                        if isOcc(v) ~= s && v ~= s
                            isOcc(v) = s;
                            break;
                        end
                    end
                    
                    n = n + 1;
                    b_surr = segY(v,:) - mean(segY(v,:));
                    
                    if std(a) < 1e-6 || std(b_surr) < 1e-6
                        continue;
                    end
                    
                    r_surr = xcorr(a, b_surr, maxlagHz, 'coeff');
                    z_surr = fisherZ(r_surr);
                    
                    % NEW: Accumulate segment-specific pseudo values
                    pseudo_signed_accum(s,:) = pseudo_signed_accum(s,:) + z_surr;
                    pseudo_abs_accum(s,:) = pseudo_abs_accum(s,:) + abs(z_surr);
                    pseudo_counts(s) = pseudo_counts(s) + 1;
                    
                    % Update meanccorrPseudoZ (accumulate Z values as in R)
                    for j = 1:lagtimes2
                        meanccorrPseudoZ(j) = meanccorrPseudoZ(j) + abs(z_surr(j));
                        meanccorrPseudoZNotAbs(j) = meanccorrPseudoZNotAbs(j) + z_surr(j);
                        meanccorrPseudo(j) = meanccorrPseudo(j) + abs(r_surr(j));
                        nPseudo(s) = nPseudo(s) + abs(r_surr(j));
                    end
                    nPseudo(s) = nPseudo(s) / (anzahlPseudosProEpoche * lagtimes2);
                end
            end
            
            % Final averages as in R
            meanccorrPseudo = meanccorrPseudo / n;
            meanccorrPseudoZ = meanccorrPseudoZ / n;
            meanccorrPseudoZNotAbs = meanccorrPseudoZNotAbs / n;
            
        else
            % METHOD 2: restrict.surrogates = FALSE (DEFAULT)
            % This is the main R algorithm
            for s = 1:numberEpochen
                a = segX(s,:) - mean(segX(s,:));
                b = segY(s,:) - mean(segY(s,:));
                
                if std(a) < 1e-6 || std(b) < 1e-6
                    continue;
                end
                
                % Real correlation for this segment (same as above)
                r_real = xcorr(a, b, maxlagHz, 'coeff');
                z_real = fisherZ(r_real);
                
                Z_real_signed(s,:) = z_real;
                Z_real_abs(s,:) = abs(z_real);
                
                % Compare with ALL other segments (including self for real)
                for v = 1:numberEpochen
                    if s == v
                        % Real synchrony (i == h in R code)
                        for j = 1:lagtimes2
                            meanccorrRealZ(j) = meanccorrRealZ(j) + abs(z_real(j));
                            meanccorrRealZNotAbs(j) = meanccorrRealZNotAbs(j) + z_real(j);
                            meanccorrReal(j) = meanccorrReal(j) + abs(r_real(j));
                            nReal(s) = nReal(s) + abs(r_real(j));
                        end
                        nReal(s) = nReal(s) / lagtimes2;
                    else
                        % Pseudo synchrony (i != h in R code)
                        b_surr = segY(v,:) - mean(segY(v,:));
                        
                        if std(a) < 1e-6 || std(b_surr) < 1e-6
                            continue;
                        end
                        
                        r_surr = xcorr(a, b_surr, maxlagHz, 'coeff');
                        z_surr = fisherZ(r_surr);
                        
                        % NEW: Accumulate segment-specific pseudo values
                        pseudo_signed_accum(s,:) = pseudo_signed_accum(s,:) + z_surr;
                        pseudo_abs_accum(s,:) = pseudo_abs_accum(s,:) + abs(z_surr);
                        pseudo_counts(s) = pseudo_counts(s) + 1;
                        
                        for j = 1:lagtimes2
                            meanccorrPseudoZ(j) = meanccorrPseudoZ(j) + abs(z_surr(j));
                            meanccorrPseudoZNotAbs(j) = meanccorrPseudoZNotAbs(j) + z_surr(j);
                            meanccorrPseudo(j) = meanccorrPseudo(j) + abs(r_surr(j));
                            nPseudo(s) = nPseudo(s) + abs(r_surr(j));
                        end
                        nPseudo(s) = nPseudo(s) / ((numberEpochen - 1) * lagtimes2);
                    end
                end
            end
            
            % Final averages as in R
            total_pseudo_comparisons = (numberEpochen - 1) * numberEpochen;
            meanccorrPseudo = meanccorrPseudo / total_pseudo_comparisons;
            meanccorrPseudoZ = meanccorrPseudoZ / total_pseudo_comparisons;
            meanccorrPseudoZNotAbs = meanccorrPseudoZNotAbs / total_pseudo_comparisons;
        end
        
        % NEW: Calculate segment-specific pseudo averages
        for s = 1:numberEpochen
            if pseudo_counts(s) > 0
                Z_pseudo_signed(s,:) = pseudo_signed_accum(s,:) / pseudo_counts(s);
                Z_pseudo_abs(s,:) = pseudo_abs_accum(s,:) / pseudo_counts(s);
            end
        end
        
        % Final averages for real (same for both methods)
        meanccorrReal = meanccorrReal / numberEpochen;
        meanccorrRealZ = meanccorrRealZ / numberEpochen;
        meanccorrRealZNotAbs = meanccorrRealZNotAbs / numberEpochen;

        % Calculate k1 and k1NotAbs (as in R)
        k1 = 0;
        k1NotAbs = 0;
        for t = 1:lagtimes2
            if meanccorrRealZ(t) > meanccorrPseudoZ(t)
                k1 = k1 + 1;
            end
            if meanccorrRealZNotAbs(t) > meanccorrPseudoZNotAbs(t)
                k1NotAbs = k1NotAbs + 1;
            end
        end

        % Calculate Effect Size (ES) as in R's as.data.frame function
        pseudo_std = std(meanccorrPseudoZ);
        if pseudo_std < 1e-10
            ES = NaN;
        else
            ES = (mean(meanccorrRealZ) - mean(meanccorrPseudoZ)) / pseudo_std;
        end
        
        pseudo_std_notAbs = std(meanccorrPseudoZNotAbs);
        if pseudo_std_notAbs < 1e-10
            ES_notAbs = NaN;
        else
            ES_notAbs = (mean(meanccorrRealZNotAbs) - mean(meanccorrPseudoZNotAbs)) / pseudo_std_notAbs;
        end
        
        % Additional ES calculations for leads (as in R)
        Z_lead1 = mean(meanccorrRealZ(1:maxlagHz));
        Z_lead2 = mean(meanccorrRealZ((maxlagHz+2):lagtimes2));
        
        if maxlagHz >= 1
            pseudo_std_lead1 = std(meanccorrPseudoZ(1:maxlagHz));
            if pseudo_std_lead1 < 1e-10
                ES_lead1 = NaN;
            else
                ES_lead1 = (mean(meanccorrRealZ(1:maxlagHz)) - mean(meanccorrPseudoZ(1:maxlagHz))) / pseudo_std_lead1;
            end
        else
            Z_lead1 = NaN;
            ES_lead1 = NaN;
        end
        
        if lagtimes2 > (maxlagHz + 1)
            pseudo_std_lead2 = std(meanccorrPseudoZ((maxlagHz+2):lagtimes2));
            if pseudo_std_lead2 < 1e-10
                ES_lead2 = NaN;
            else
                ES_lead2 = (mean(meanccorrRealZ((maxlagHz+2):lagtimes2)) - mean(meanccorrPseudoZ((maxlagHz+2):lagtimes2))) / pseudo_std_lead2;
            end
        else
            Z_lead2 = NaN;
            ES_lead2 = NaN;
        end

        % For non-absolute values (in-phase/anti-phase)
        positive_idx = meanccorrRealZNotAbs >= 0;
        negative_idx = meanccorrRealZNotAbs < 0;
        meanZ_inphase = mean(meanccorrRealZNotAbs(positive_idx));
        meanZ_antiphase = mean(meanccorrRealZNotAbs(negative_idx));
        count_inphase = sum(positive_idx);
        count_antiphase = sum(negative_idx);

        % Store summary statistics (matches R's as.data.frame output)
        summary_row = table();
        summary_row.Var1 = ids(pairs(p,1));
        summary_row.Var2 = ids(pairs(p,2));
        summary_row.n_data = size_val;
        summary_row.Z = mean(meanccorrRealZ);
        summary_row.Z_Pseudo = mean(meanccorrPseudoZ);
        summary_row.SD_Z = std(meanccorrRealZ);
        summary_row.SD_Z_Pseudo = std(meanccorrPseudoZ);
        summary_row.n_lags = lagtimes2;
        summary_row.Perc_gt_Pseudo = 100 * k1 / lagtimes2;
        summary_row.n_Segments = numberEpochen;
        summary_row.ES = ES;
        summary_row.Z_lead1 = Z_lead1;
        summary_row.Z_lead2 = Z_lead2;
        summary_row.ES_lead1 = ES_lead1;
        summary_row.ES_lead2 = ES_lead2;
        summary_row.meanZ_inphase = meanZ_inphase;
        summary_row.meanZ_antiphase = meanZ_antiphase;
        summary_row.count_inphase = count_inphase;
        summary_row.count_antiphase = count_antiphase;
        summary_row.Z_noAbs = mean(meanccorrRealZNotAbs);
        summary_row.Z_Pseudo_noAbs = mean(meanccorrPseudoZNotAbs);
        summary_row.Perc_gt_Pseudo_noAbs = 100 * k1NotAbs / lagtimes2;
        summary_row.ES_noAbs = ES_notAbs;
        
        summary_stats = [summary_stats; summary_row];
        
        % Store ES values in ES_table
        es_row = table();
        es_row.Dyad1 = ids(pairs(p,1));
        es_row.Dyad2 = ids(pairs(p,2));
        es_row.ES = ES;
        es_row.ES_pseudo = mean(meanccorrPseudoZ);
        es_row.ES_noAbs = ES_notAbs;
        es_row.ES_lead1 = ES_lead1;
        es_row.ES_lead2 = ES_lead2;
        
        ES_table = [ES_table; es_row];

        % Create detailed output table (per segment, per lag)
        num_entries = numberEpochen * lagtimes2;
        
        dyad1_col = repmat(ids(pairs(p,1)), num_entries, 1);
        dyad2_col = repmat(ids(pairs(p,2)), num_entries, 1);
        segment_col = repelem((1:numberEpochen)', lagtimes2);
        lag_col = repmat(lags(:), numberEpochen, 1);
        
        % Reshape matrices - NOW USING SEGMENT-SPECIFIC VALUES
        z_col = reshape(Z_real_signed', [], 1);
        z_abs_col = reshape(Z_real_abs', [], 1);
        z_pseudo_col = reshape(Z_pseudo_signed', [], 1);
        z_pseudo_abs_col = reshape(Z_pseudo_abs', [], 1);

        temp = table(dyad1_col, dyad2_col, segment_col, lag_col, ...
                     z_col, z_abs_col, z_pseudo_col, z_pseudo_abs_col, ...
                    'VariableNames', {'Dyad1','Dyad2','Segment','Lag',...
                                      'Z','Z_abs','Z_pseudo','Z_pseudo_abs'});

        tbl = [tbl; temp];
    end
end

function z = fisherZ(r)
    % Fisher Z transform (matches R's 0.5*log((1+r)/(1-r)))
    r = max(min(r, 0.999), -0.999);
    z = 0.5 * log((1 + r) ./ (1 - r));
end

function output = mv_susy_core(X, cfg)
% MV_SUSY_CORE - Compute multivariate synchrony with surrogate testing

% === Parameters ===
method = getfieldwithdefault(cfg, 'method', 'lambda_max');  % or 'omega'
Hz = getfieldwithdefault(cfg, 'fs', 1);
seg_len = getfieldwithdefault(cfg, 'segment', 5);  % in seconds
n_pseudo = getfieldwithdefault(cfg, 'total_surrogates', 500);
rng_seed = getfieldwithdefault(cfg, 'seed', 1);

if ~ismatrix(X) || size(X,2) < 2
    error('Input X must be T×N with at least 2 participants.');
end

if ~ismember(method, {'lambda_max', 'omega'})
    error('Unknown mvSUSY method: "%s". Use "lambda_max" or "omega".', method);
end

if ~isempty(rng_seed)
    rng(rng_seed);
end

% === Segment the data ===
T = size(X,1);
N = size(X,2);
segment_Hz = seg_len * Hz;
nSeg = floor(T / segment_Hz);

if nSeg < 3
    error('Too few segments (%d). Reduce segment length or increase data.', nSeg);
end

X = X(1:nSeg*segment_Hz, :);  % truncate
X = reshape(X, [segment_Hz, nSeg, N]);  % time × segment × subject

% === Define synchrony metric ===
switch lower(method)
    case 'lambda_max'
        metric_fun = @(M) max(eig(corrcoef(M, 'Rows', 'pairwise')));
    case 'omega'
        metric_fun = @(M) 1 - det(abs(cov(M))) / prod(diag(cov(M)));
end

% === Compute real synchrony scores ===
sync_real = nan(nSeg,1);
for s = 1:nSeg
    segMat = squeeze(X(:,s,:));  % T × N
    sync_real(s) = metric_fun(segMat);
end

% === Generate pseudo segments ===
pseudo_matrix = nan(n_pseudo, nSeg);
for p = 1:n_pseudo
    permuted = X(:,randperm(nSeg),:);  % scramble segment order
    for s = 1:nSeg
        pseudo_matrix(p,s) = metric_fun(squeeze(permuted(:,s,:)));
    end
end

% === Compile results ===
output = struct();
output.method = method;
output.real = sync_real;
output.pseudo = pseudo_matrix;
output.mean_real = mean(sync_real);
output.mean_pseudo = mean(pseudo_matrix(:));
output.effect_size = (output.mean_real - output.mean_pseudo) / std(pseudo_matrix(:));
[~, output.p_value] = ttest2(sync_real, pseudo_matrix(:));
output.segmentHz = segment_Hz;
output.nSeg = nSeg;
output.nPseudo = n_pseudo;
end

function val = getfieldwithdefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

function T = mv_susy_summary_table(res)
% MV_SUSY_SUMMARY_TABLE: Format mvSUSY result into table for t_susy

T = table();

T.method        = string(res.method);
T.n_segments    = res.nSeg;
T.segment_Hz    = res.segmentHz;
T.n_pseudo      = res.nPseudo;

T.mean_real     = res.mean_real;
T.mean_pseudo   = res.mean_pseudo;
T.effect_size   = res.effect_size;

T.p_value       = res.p_value;  % Changed from t_pval to p_value

end