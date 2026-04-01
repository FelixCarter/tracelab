function summary_data = t_movement_summary(dat, cfg)
% T_MOVEMENT_SUMMARY Summarizes head movement data across participants and conditions

    % Set default configuration
    if nargin < 2
        cfg = struct();
    end
    
    if ~isfield(cfg, 'avg_ppt'), cfg.avg_ppt = true; end
    if ~isfield(cfg, 'fs'), cfg.fs = 1; end
    if ~isfield(cfg, 'bin'), cfg.bin = 1; end
    
    % Get condition names
    condition_names = fieldnames(dat);
    num_conditions = length(condition_names);
    
    % Initialize output structure
    summary_data = struct();
    summary_data.cfg = cfg;
    summary_data.condition_names = condition_names;
    
    % Process each condition
    for cond_idx = 1:num_conditions
        cond_name = condition_names{cond_idx};
        participants = dat.(cond_name);
        num_participants = length(participants);
        
        % Initialize storage for this condition
        condition_data = struct();
        condition_data.participant_tables = cell(1, num_participants);
        condition_data.has_valid_data = false(1, num_participants);
        
        % Process each participant in this condition
        for ppt_idx = 1:num_participants
            participant_data = participants{ppt_idx};
            
            if isfield(participant_data, 'HeadMovement') && isfield(participant_data, 'HeadMovementTime')
                head_movement = participant_data.HeadMovement;
                head_movement_time = participant_data.HeadMovementTime;
                
                % Create binned data
                [time_bins, movement_bins] = bin_movement_data(head_movement, head_movement_time, cfg);
                
                if ~isempty(time_bins)
                    % Store participant data as table
                    condition_data.participant_tables{ppt_idx} = table(time_bins, movement_bins, ...
                        'VariableNames', {'Time', 'Movement'});
                    condition_data.has_valid_data(ppt_idx) = true;
                else
                    condition_data.participant_tables{ppt_idx} = table();
                    condition_data.has_valid_data(ppt_idx) = false;
                end
            else
                warning('Participant %d in condition %s missing HeadMovement or HeadMovementTime fields', ...
                    ppt_idx, cond_name);
                condition_data.participant_tables{ppt_idx} = table();
                condition_data.has_valid_data(ppt_idx) = false;
            end
        end
        
        % Store condition data
        summary_data.(cond_name) = condition_data;
    end
    
    % If averaging across participants is requested
    if cfg.avg_ppt
        summary_data = create_averaged_data(summary_data);
    end
end

function [time_bins, movement_bins] = bin_movement_data(head_movement, head_movement_time, cfg)
% Bin movement data into time segments
    
    if isempty(head_movement) || isempty(head_movement_time)
        time_bins = [];
        movement_bins = [];
        return;
    end
    
    % Ensure vectors are columns
    head_movement = head_movement(:);
    head_movement_time = head_movement_time(:);
    
    % Calculate bin edges
    min_time = min(head_movement_time);
    max_time = max(head_movement_time);
    bin_edges = min_time:cfg.bin:max_time;
    
    if isempty(bin_edges) || length(bin_edges) < 2
        time_bins = mean(head_movement_time);
        movement_bins = mean(head_movement);
        return;
    end
    
    % Add final edge if needed
    if bin_edges(end) < max_time
        bin_edges(end+1) = bin_edges(end) + cfg.bin;
    end
    
    num_bins = length(bin_edges) - 1;
    time_bins = zeros(num_bins, 1);
    movement_bins = zeros(num_bins, 1);
    
    % Calculate mean movement for each bin
    for bin_idx = 1:num_bins
        bin_start = bin_edges(bin_idx);
        bin_end = bin_edges(bin_idx + 1);
        
        % Find samples in this bin
        in_bin = head_movement_time >= bin_start & head_movement_time < bin_end;
        
        if any(in_bin)
            time_bins(bin_idx) = (bin_start + bin_end) / 2;
            movement_bins(bin_idx) = mean(head_movement(in_bin), 'omitnan');
        else
            time_bins(bin_idx) = (bin_start + bin_end) / 2;
            movement_bins(bin_idx) = NaN;
        end
    end
    
    % Remove bins with no data
    valid_bins = ~isnan(movement_bins);
    time_bins = time_bins(valid_bins);
    movement_bins = movement_bins(valid_bins);
end

function summary_data = create_averaged_data(summary_data)
% Create averaged data across participants for each condition
    
    condition_names = summary_data.condition_names;
    
    for cond_idx = 1:length(condition_names)
        cond_name = condition_names{cond_idx};
        condition_data = summary_data.(cond_name);
        participant_tables = condition_data.participant_tables;
        has_valid_data = condition_data.has_valid_data;
        
        if ~any(has_valid_data)
            condition_data.averaged_table = table();
            summary_data.(cond_name) = condition_data;
            continue;
        end
        
        % Use only participants with valid data
        valid_tables = participant_tables(has_valid_data);
        
        % Find common time points across all participants
        all_times = [];
        for ppt_idx = 1:length(valid_tables)
            if ~isempty(valid_tables{ppt_idx})
                all_times = union(all_times, valid_tables{ppt_idx}.Time);
            end
        end
        
        if isempty(all_times)
            condition_data.averaged_table = table();
            summary_data.(cond_name) = condition_data;
            continue;
        end
        
        % Initialize arrays for averaging
        movement_matrix = NaN(length(all_times), length(valid_tables));
        
        % Collect movement data for each time point
        for ppt_idx = 1:length(valid_tables)
            ppt_table = valid_tables{ppt_idx};
            for time_idx = 1:length(all_times)
                match_idx = find(abs(ppt_table.Time - all_times(time_idx)) < 1e-10, 1); % tolerance against floating point rounding errors
                if ~isempty(match_idx)
                    movement_matrix(time_idx, ppt_idx) = ppt_table.Movement(match_idx);
                end
            end
        end
        
        % Calculate mean and standard error
        avg_movement = mean(movement_matrix, 2, 'omitnan');
        sem_movement = std(movement_matrix, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(movement_matrix), 2));
        
        % Store averaged data
        condition_data.averaged_table = table(all_times, avg_movement, sem_movement, ...
            'VariableNames', {'Time', 'MeanMovement', 'SEM'});
        summary_data.(cond_name) = condition_data;
    end
end