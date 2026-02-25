function preproc_output = t_preproc(cfg)

%--------------------------------------------------------------------------
% Configuration Parameters for t_preproc(cfg)
%
% Required:
%   cfg.dir                   → Path to input .csv directory

% Optional (defaults listed):
%   cfg.groupvar              → Grouping variable name (default: "")
%   cfg.unify                 → How to unify durations across groups:
%                                 'modal' | 'median' | 'mean' | 'shortest' | 'longest'  (default: 'modal')
%   cfg.round                 → round 'up' or 'down' when unifying data lengths?
%   cfg.fill                  → Gap-filling method for missing frames:
%                                 'interp' only currently (default: 'interp')
%   cfg.target_fs             → Target sampling frequency in Hz (default: 5)
%   cfg.bin                   → Binning toggle for landmark data (default: 0 = off)
%   cfg.nose_end              → Landmark index where nose ends (default: 34)
%   cfg.motion_window         → Smoothing window for head motion in seconds (default: 1.0)
%   cfg.emo_bin               → Smoothing window for emotions in seconds (default: 0.2)
%   cfg.min_valid_rows        → Minimum valid rows per group (default: 5)
%   cfg.min_samples_per_second→ Minimum samples/sec per valid time window (default: 2)
%   cfg.max_bad_seconds_ratio → Max proportion of under-sampled seconds before skipping (default: 0.5)
%   cfg.remove_motionless     → Whether to remove segments with zero motion variance (default: true)
%   cfg.participant_id_column → Column name for participant ID (default: 'participant')
%--------------------------------------------------------------------------

% --- Set Defaults ---
if ~isfield(cfg, 'dir'),        error('cfg.dir must be provided'); end
if ~isfield(cfg, 'groupvar'),   cfg.groupvar = ""; end
if ~isfield(cfg, 'unify'),      cfg.unify = 'modal'; end
if ~isfield(cfg, 'fill'),       cfg.fill = 'interp'; end
if ~isfield(cfg, 'target_fs'),  cfg.target_fs = 5; end
if ~isfield(cfg, 'bin'),        cfg.bin = 0; end
if ~isfield(cfg, 'nose_end'),   cfg.nose_end = 34; end
if ~isfield(cfg, 'motion_window'), cfg.motion_window = 1.0; end
if ~isfield(cfg, 'emo_bin') || isempty(cfg.emo_bin), cfg.emo_bin = 0.2; end
if ~isfield(cfg, 'min_valid_rows'), cfg.min_valid_rows = 5; end
if ~isfield(cfg, 'min_samples_per_second'), cfg.min_samples_per_second = 2; end
if ~isfield(cfg, 'max_bad_seconds_ratio'), cfg.max_bad_seconds_ratio = 0.5; end
if ~isfield(cfg, 'remove_motionless'), cfg.remove_motionless = true; end
if ~isfield(cfg, 'participant_id_column'), cfg.participant_id_column = 'participant'; end
if ~isfield(cfg, 'round'), cfg.round = 'down'; end

files = dir(fullfile(cfg.dir, '*.csv'));
groupDurations = containers.Map();
participantStore = struct();

% --- First pass: gather durations per group ---
for f = 1:length(files)
    dat = readtable(fullfile(cfg.dir, files(f).name));
    requiredCols = {'time','landmarks', cfg.participant_id_column};
    if ~all(ismember(requiredCols, dat.Properties.VariableNames)), continue; end

    groupvar = cfg.groupvar;
    if groupvar == "" || ~ismember(groupvar, dat.Properties.VariableNames)
        dat.groupDummy = repmat("allData", height(dat), 1); groupvar = 'groupDummy';
    end
    dat.(groupvar) = string(dat.(groupvar));
    groups = unique(dat.(groupvar));
    groups = groups(groups ~= "");

    for i = 1:length(groups)
        rows = strcmp(dat.(groupvar), groups(i));
        tVals = dat.time(rows);
        if strcmpi(cfg.round, 'up')
            dur = ceil(max(tVals) - min(tVals));
        else
            dur = floor(max(tVals) - min(tVals));
        end
        k = matlab.lang.makeValidName(groups(i));
        if ~isKey(groupDurations, k), groupDurations(k) = []; end
        groupDurations(k) = [groupDurations(k), dur];
    end
end

% --- Compute unified durations ---
unifiedDur = struct();
keysList = keys(groupDurations);
for i = 1:length(keysList)
    k = keysList{i}; durs = groupDurations(k);
    switch cfg.unify
        case 'modal',    unifiedDur.(k) = mode(durs);
        case 'median',   unifiedDur.(k) = median(durs);
        case 'mean',     unifiedDur.(k) = mean(durs);
        case 'shortest', unifiedDur.(k) = min(durs);
        case 'longest',  unifiedDur.(k) = max(durs);
    end
end

% --- Second pass: process each file ---
for f = 1:length(files)
    fname = files(f).name;
    msgPrefix = sprintf('[%d/%d]', f, length(files));
    dat = readtable(fullfile(cfg.dir, fname));

    required = {'time','landmarks', cfg.participant_id_column};
    hasGroup = cfg.groupvar == "" || ismember(cfg.groupvar, dat.Properties.VariableNames);
    if ~all(ismember(required, dat.Properties.VariableNames)) || ~hasGroup
        fprintf('%s Skipped: %s (missing required columns)\n', msgPrefix, fname);
        continue;
    end

    groupvar = cfg.groupvar;
    if groupvar == "" || ~ismember(groupvar, dat.Properties.VariableNames)
        dat.groupDummy = repmat("allData", height(dat), 1); groupvar = 'groupDummy';
    end
    dat.(groupvar) = string(dat.(groupvar));
    groups = unique(dat.(groupvar));
    groups = groups(groups ~= "");
    processedSomething = false;

    for c = 1:length(groups)
        try
            rows = strcmp(dat.(groupvar), groups(c));
            t = dat.time(rows);
            if length(t) < 2, continue; end

            % Per-second sample filtering
            timeRaw = dat.time(rows);
            secondLabels = floor(timeRaw);
            [secs, ~, secIdx] = unique(secondLabels);
            countsPerSec = accumarray(secIdx, 1);
            badSecs = countsPerSec < cfg.min_samples_per_second;
            badRatio = sum(badSecs) / numel(secs);

            if badRatio >= cfg.max_bad_seconds_ratio
                fprintf('%s Skipped: %s → Group: %s (%.1f%% bad seconds)\n', ...
                    msgPrefix, fname, groups(c), 100 * badRatio);
                continue;
            end

            key = matlab.lang.makeValidName(groups(c));
            if ~isfield(unifiedDur, key), continue; end

            targetT = (0:1/cfg.target_fs:unifiedDur.(key))';
            [xyStacked, xyTable, valid] = process_landmarks(dat, rows, targetT, cfg);
            if isempty(xyStacked), continue; end

            out = struct();
            [out.HeadMovement, out.HeadMovementTime] = compute_head_motion(xyStacked, targetT, cfg);

            % Optional emotion data
            emoVars = {'happyDetection','sadDetection','neutralDetection','angryDetection',...
                       'disgustedDetection','surprisedDetection','fearfulDetection'};
            if all(ismember(emoVars, dat.Properties.VariableNames))
                emo = dat{rows, emoVars};
                emo = emo(valid,:);
                minLen = min(length(t), size(emo,1));
                t = t(1:minLen); emo = emo(1:minLen,:);
                validEmo = ~isnan(t) & all(~isnan(emo), 2);
                tEmo = t(validEmo); emo = emo(validEmo,:);
                if numel(tEmo) >= 2
                    emoTable = process_emotions(tEmo, emo, targetT, cfg);
                else
                    emoTable = array2table(NaN(length(targetT), length(emoVars)), 'VariableNames', emoVars);
                    emoTable = [table(targetT, 'VariableNames', {'time'}), emoTable];
                end
                out.Emotions = emoTable;
            end

            % Use user-defined participant ID column
            if ~ismember(cfg.participant_id_column, dat.Properties.VariableNames)
                error('Participant ID column "%s" not found in %s', cfg.participant_id_column, fname);
            end
            out.ParticipantID = string(dat.(cfg.participant_id_column)(find(rows, 1, 'first')));
            out.Time = targetT;
            out.XY = xyTable;
            % Store the original filename in Settings
            settingsWithFilename = cfg;
            settingsWithFilename.source_filename = fname;  % Add filename to settings
            out.Settings = settingsWithFilename;

            if ~isfield(participantStore, key)
                participantStore.(key) = {};
            end
            participantStore.(key){end+1} = out;

            fprintf('%s Processed: %s → Group: %s\n', msgPrefix, fname, groups(c));
            save('partial_output.mat', 'participantStore');
            processedSomething = true;

        catch ME
            fprintf('%s Error: %s → Group: %s\n     ↪ %s\n', msgPrefix, fname, groups(c), ME.message);
        end
    end

    if ~processedSomething
        fprintf('%s Skipped: %s (no valid groups or landmark data)\n', msgPrefix, fname);
    end
end

% Remove motionless
if cfg.remove_motionless
    clips = fieldnames(participantStore);
    for i = 1:numel(clips)
        clip = clips{i};
        entries = participantStore.(clip);
        for j = numel(entries):-1:1
            p = entries{j};
            if ~isfield(p, 'HeadMovement') || isempty(p.HeadMovement) ...
               || std(p.HeadMovement, 'omitnan') == 0
                fprintf('Dropped motionless participant: %s (clip: %s)\n', ...
                    p.ParticipantID, clip);
                entries(j) = [];
            end
        end
        participantStore.(clip) = entries;
    end
end

preproc_output = participantStore;
end

function [xyStacked, xyTable, valid] = process_landmarks(dat, rows, targetT, cfg)
    if ~isfield(cfg, 'min_valid_rows'), cfg.min_valid_rows = 5; end
    if ~isfield(cfg, 'min_frames_per_sec'), cfg.min_frames_per_sec = 2; end
    if ~isfield(cfg, 'min_sec_coverage'), cfg.min_sec_coverage = 0.5; end

    subset = dat(rows, :);
    landmarkStrs = subset.landmarks;
    isValid = false(height(subset), 1);
    coords = cell(height(subset), 1);

    for i = 1:height(subset)
        try
            val = jsondecode(landmarkStrs{i});
            if isnumeric(val) && size(val,1) == 68 && size(val,2) == 2
                isValid(i) = true;
                coords{i} = val;
            end
        catch
            % ignore malformed rows
        end
    end

    % Trim to valid
    cleanT = subset.time(isValid);
    coords = coords(isValid);
    valid = isValid;

    if numel(cleanT) < cfg.min_valid_rows
        xyStacked = []; xyTable = []; return;
    end

    % Stack as [N × 136]
    xyStacked = cell2mat(cellfun(@(xy) reshape(xy', 1, []), coords, 'UniformOutput', false));

    % Deduplicate timestamps
    [uniqueT, ~, idxMap] = unique(cleanT);
    if length(uniqueT) < length(cleanT)
        xyStacked = grpstats(xyStacked, idxMap, 'mean');
        cleanT = uniqueT;
    end

    roundedSec = floor(cleanT);
    secCounts = accumarray(roundedSec - min(roundedSec) + 1, 1);
    validSecs = sum(secCounts >= cfg.min_frames_per_sec);
    totalSecs = max(roundedSec) - min(roundedSec) + 1;
    coverage = validSecs / totalSecs;
    fprintf('ℹ️ %d of %d seconds had ≥%d valid frames (%.1f%% coverage)\n', ...
    validSecs, totalSecs, cfg.min_frames_per_sec, 100 * coverage);

    if coverage < cfg.min_sec_coverage
        fprintf('Skipped due to sparse temporal coverage: %.1f%% of seconds had data\n', 100 * coverage);
        xyStacked = []; xyTable = []; valid = [];
        return;
    end


    try
        interpXY = interp1(cleanT, xyStacked, targetT, 'linear', 'extrap');
    catch
        warning('interp1 failed for current segment — skipping');
        xyStacked = []; xyTable = []; return;
    end

    % Labels
    xyLabels = reshape([ ...
        arrayfun(@(i) sprintf('X_%02d', i), 1:68, 'UniformOutput', false); ...
        arrayfun(@(i) sprintf('Y_%02d', i), 1:68, 'UniformOutput', false) ...
    ], 1, []);

    xyTable = array2table(interpXY, 'VariableNames', xyLabels);
    xyTable.time = targetT;
    xyTable = movevars(xyTable, 'time', 'Before', 1);
end

function [motion_cm_per_s, T1Hz] = compute_head_motion(xyStacked, targetT, cfg)
    % Downsample to 1 Hz based on rounded time
    roundedTime = round(targetT);
    [T1Hz, uniqIdx] = unique(roundedTime, 'stable');
    uniqIdx = min(uniqIdx, size(xyStacked, 1));  % cap to size
    xy1Hz = xyStacked(uniqIdx, :);

    % Eye landmark indices
    leftEyeIdx  = sort([2*(37:42)-1, 2*(37:42)]);
    rightEyeIdx = sort([2*(43:48)-1, 2*(43:48)]);

    leftX  = mean(xy1Hz(:, leftEyeIdx(1:2:end)), 2);
    leftY  = mean(xy1Hz(:, leftEyeIdx(2:2:end)), 2);
    rightX = mean(xy1Hz(:, rightEyeIdx(1:2:end)), 2);
    rightY = mean(xy1Hz(:, rightEyeIdx(2:2:end)), 2);

    IOD_px = sqrt((rightX - leftX).^2 + (rightY - leftY).^2);
    IOD_px(IOD_px == 0) = NaN;

    % Nose tip motion
    noseIdx = cfg.nose_end;
    xCol = 2 * (noseIdx - 1) + 1;
    yCol = xCol + 1;
    noseX = xy1Hz(:, xCol);
    noseY = xy1Hz(:, yCol);

    dx = [0; diff(noseX)];
    dy = [0; diff(noseY)];
    motion_px = sqrt(dx.^2 + dy.^2);

    % Convert to cm/s
    IOD_cm = 6.3;
    motion_cm_per_s = motion_px .* (IOD_cm ./ IOD_px);
    motion_cm_per_s(1) = 0;
    motion_cm_per_s(motion_cm_per_s > 20) = NaN;
end

function emoTable = process_emotions(tRaw, emoRaw, targetT, cfg)
% process_emotions – Bin and impute emotion scores with LOCF and NOCB
%
% INPUTS:
%   tRaw     = raw timestamps associated with emotion samples
%   emoRaw   = raw emotion scores (must match tRaw in length)
%   targetT  = unified time vector from main pipeline
%   cfg      = configuration structure (expects field 'emo_bin')
%
% OUTPUT:
%   emoTable = table with 'time' column followed by emotion columns

    % Define bins (based on cfg.emo_bin)
    bins = 0:cfg.emo_bin:max(targetT) + cfg.emo_bin;
    if numel(bins) < 2
        bins = [0, max(targetT) + cfg.emo_bin];
    end
    emoBinC = bins(1:end-1);  % Use lower bin edges as timestamps

    % Initialize binned matrix
    emoBinned = NaN(length(emoBinC), size(emoRaw,2));

    % Bin means (omit NaNs)
    for b = 1:length(emoBinC)
        idx = tRaw >= bins(b) & tRaw < bins(b+1);
        if any(idx)
            emoBinned(b,:) = mean(emoRaw(idx,:), 1, 'omitnan');
        end
    end

    % LOCF – Last Observation Carried Forward
    for j = 1:size(emoBinned,2)
        for i = 2:size(emoBinned,1)
            if isnan(emoBinned(i,j))
                emoBinned(i,j) = emoBinned(i-1,j);
            end
        end
    end

    % NOCB – Fill leading NaNs
    for j = 1:size(emoBinned,2)
        if isnan(emoBinned(1,j))
            firstValid = find(~isnan(emoBinned(:,j)), 1, 'first');
            if ~isempty(firstValid)
                emoBinned(1:firstValid-1,j) = emoBinned(firstValid,j);
            end
        end
    end

    % Output table with emotion scores and timestamps
    emoVars = {'happyDetection','sadDetection','neutralDetection','angryDetection',...
               'disgustedDetection','surprisedDetection','fearfulDetection'};

    emoTable = array2table(emoBinned, 'VariableNames', emoVars);
    emoTable = [table(emoBinC', 'VariableNames', {'time'}), emoTable];
end