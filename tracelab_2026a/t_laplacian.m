function dat_laplacian = t_laplacian(dat, meanFace, cfg)
    % Laplacian filtering: mean-center landmarks within facial regions

    % Default region map if not provided
    default_map = containers.Map;
    default_map('group_jawline')   = 1:17;
    default_map('group_right_eye') = 37:42;
    default_map('group_left_eye')  = 43:48;
    default_map('group_nose')      = 27:36;
    default_map('group_mouth')     = 49:68;

    % Use cfg.region_map if given, otherwise default
    if nargin < 3 || ~isfield(cfg, 'region_map')
        region_map = default_map;
    else
        region_map = cfg.region_map;
    end

    dat_laplacian = dat;
    conditions = fieldnames(dat);

    for c = 1:numel(conditions)
        condName = conditions{c};
        dataCells = dat.(condName);
        nSubj = numel(dataCells);

        for s = 1:nSubj
            subjData = dataCells{s};
            XYtable = subjData.XY;
            T = size(XYtable, 1);
            rawXY = table2array(XYtable(:,2:end));         % [T × D]
            laplaceFrameMatrix = rawXY;                    % placeholder

            for region = region_map.keys
                idxs = region_map(region{1});
                xyIdxs = reshape([2*idxs - 1; 2*idxs], 1, []);  % X and Y column indices

                regionXY = rawXY(:, xyIdxs);                    % [T × 2N]
                coords = reshape(regionXY', 2, [], T);          % [2 × landmarks × T]
                meanCoords = mean(coords, 2);                   % [2 × 1 × T]
                lap = coords - meanCoords;                      % [2 × landmarks × T]
                lap = reshape(lap, [], T)';                     % [T × 2N]

                laplaceFrameMatrix(:, xyIdxs) = lap;
            end

            % Preserve original structure
            timeVec = XYtable{:,1};
            varNames = XYtable.Properties.VariableNames;
            lapTable = array2table([timeVec, laplaceFrameMatrix], 'VariableNames', varNames);

            subjData.XY = lapTable;
            dat_laplacian.(condName){s} = subjData;
        end
    end
end
