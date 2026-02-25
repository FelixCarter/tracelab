function summary_out = t_susy_summarise(susy_out, cfg)
% t_susy_summarise - Hierarchical averaging of SUSY output
% takes output from t_susy as input

    if ~isfield(cfg, 'avg_segment'), cfg.avg_segment = true; end
    if ~isfield(cfg, 'avg_dyad'),    cfg.avg_dyad    = true; end
    if ~isfield(cfg, 'avg_ppt'),     cfg.avg_ppt     = false; end

    cfg.avg_segment = is_truthy(cfg.avg_segment);
    cfg.avg_dyad    = is_truthy(cfg.avg_dyad);
    cfg.avg_ppt     = is_truthy(cfg.avg_ppt);

    summary_out = struct();
    groups = fieldnames(susy_out);

    for g = 1:numel(groups)
        group = groups{g};
        T = susy_out.(group);
        if isempty(T) || height(T) == 0
            summary_out.(group) = table(); continue;
        end

        % average across lags first (matches t_plot_susy)
        if ismember('Lag', T.Properties.VariableNames)
            vars = get_value_vars(T);
            T = groupsummary(T, {'Dyad1','Dyad2','Segment'}, ...
                @(x) mean(x, 'omitnan'), vars);
            T.Properties.VariableNames = clean_colnames(T.Properties.VariableNames);
        end

        % Step 2: Average across segments per dyad
        if cfg.avg_dyad && all(ismember({'Dyad1','Dyad2'}, T.Properties.VariableNames))
            T = removevars(T, intersect({'Segment','Lag','GroupCount'}, T.Properties.VariableNames));
            vars = get_value_vars(T);
            T = groupsummary(T, {'Dyad1','Dyad2'}, ...
                @(x) mean(x, 'omitnan'), vars);
            T.Properties.VariableNames = clean_colnames(T.Properties.VariableNames);
        end

        % Step 3: Average across dyads per participant
        if cfg.avg_ppt
            if all(ismember({'Dyad1','Dyad2'}, T.Properties.VariableNames))
                T1 = T; T1.Participant = string(T1.Dyad1);
                T2 = T; T2.Participant = string(T2.Dyad2);
                T = [T1; T2];
            end
            T = removevars(T, intersect({'Dyad1','Dyad2','Segment','Lag','GroupCount'}, T.Properties.VariableNames));
            vars = get_value_vars(T);
            T = groupsummary(T, 'Participant', ...
                @(x) mean(x, 'omitnan'), vars);
            T.Properties.VariableNames = clean_colnames(T.Properties.VariableNames);
        end

        % Final cleanup
        dropCols = intersect({'GroupCount'}, T.Properties.VariableNames);
        T = removevars(T, dropCols);

        summary_out.(group) = T;
    end
end

function out = is_truthy(val)
    if islogical(val)
        out = val;
    elseif isnumeric(val)
        out = val ~= 0;
    elseif ischar(val) || isstring(val)
        out = any(strcmpi(string(val), {'true','yes','1'}));
    else
        out = false;
    end
end

function vars = get_value_vars(T)
    candidates = {'Z','Z_abs','Z_pseudo','Z_pseudo_abs'};
    vars = T.Properties.VariableNames;
    vars = vars(contains(vars, candidates));
end

function names_out = clean_colnames(names_in)
    names_out = regexprep(names_in, '^(fun\d+_)+', '');
    names_out = regexprep(names_out, '^(mean_)+', '');
end