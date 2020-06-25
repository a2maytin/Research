function obj = fitMSDandrew(obj)
%%Fit all MSD curves by a linear function. Modified version of fitMSD.
%
% obj = obj.fitMSDandrew fits all MSD curves by a straight line
%                      y = a * x + b.
% The fit is therefore rigorously valid only for purely
% diffusive behavior.
%
% Results are stored in the 'fit' field of the returned
% object. It is a structure with 2 fields:
% - a: all the values of the slope of the linear fit.
% - b: all the values for the intersect of the linear fit.
% - r2fit: the adjusted R2 value as a indicator of the goodness
% of the fit.

% Fits first 2 points in MSD vs. tau, plus the origin. 

n_spots = numel(obj.msd);

fprintf('Fitting %d curves of MSD = f(t), taking only the first 2 points of each curve and thru the origin... ', n_spots)

a = NaN(n_spots, 1);
b = NaN(n_spots, 1);
r2fit = NaN(n_spots, 1);
ft = fittype('poly1');

fprintf('%5d/%5d', 0, n_spots);
for i_spot = 1 : n_spots %loop through every track
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d/%5d', i_spot, n_spots);
    
    msd_spot = obj.msd{i_spot};
    
    t = msd_spot(:,1); %times
    y = msd_spot(:,2); %msd
    w = msd_spot(:,4); %weights
    
    % Clip data, take first 3 points
    t_limit = 1 : 3;

    t = t(t_limit);
    y = y(t_limit);
    w = w(t_limit);
    
    % Thrash bad data
    nonnan = ~isnan(y);
    x = t(nonnan);
    y = y(nonnan);
    w = w(nonnan);
    
    if numel(y) < 2
        continue
    end
    
    fit = fitlm(x,y,'Intercept',false);
    
    a(i_spot) = fit.Coefficients{1,1};
    r2fit(i_spot) = fit.Rsquared.Ordinary;
    
end
fprintf('\b\b\b\b\b\b\b\b\b\bDone.\n')

obj.lfit = struct(...
    'a', a, ...
    'r2fit', r2fit);

end