function msmsd = getEATAMSD(obj, indices)
%%Compute the unweighted mean of all MSD curves.
%
% Results are returned as a N x 4 double array, and ordered as
% following: [ dT M ] with:
% - dT the delay vector
% - M the mean of MSD for each delay
%
% msd = obj.getEAMeanMSD(indices) only takes into account the MSD
% curves with the specified indices.

if ~obj.msd_valid
    obj = obj.computeMSD(indices);
end

if nargin < 2 || isempty(indices)
    indices = 1 : numel(obj.msd);
end

n_tracks = numel(indices);

% First, collect all possible delays
all_delays = cell(n_tracks, 1);
for i = 1 : n_tracks
    index = indices(i);
    
    if isempty( obj.msd{index} )
        continue
    end
    all_delays{i} = obj.msd{index}(:,1);
end
delays = unique( vertcat( all_delays{:} ) );
n_delays = numel(delays);

% Collect
sum_weight = zeros(n_delays, 1); %all weighted equally
sum_mean   = zeros(n_delays, 1);

% 1st pass
for i = 1 : n_tracks
    
    index = indices(i);
    if isempty( obj.msd{index} )
        continue
    end
    
    t = obj.msd{index}(:,1);
    m = obj.msd{index}(:,2);
    
    % Do not tak NaNs
    valid = ~isnan(m);
    t = t(valid);
    m = m(valid);
    
    % Find common indices
    [~, index_in_all_delays, ~] = intersect(delays, t);
    
    % Accumulate
    sum_weight(index_in_all_delays)  = sum_weight(index_in_all_delays) + 1; 
    sum_mean(index_in_all_delays)    = sum_mean(index_in_all_delays)  + m;
end

% Compute weighted mean
mmean = sum_mean ./ sum_weight;

% Output [ T mean ]
msmsd = [ delays mmean ];

end