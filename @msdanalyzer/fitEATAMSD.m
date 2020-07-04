function  varargout = fitMeanMSD(obj, clip_factor)
%%FITMEANMSD Fit the weighted averaged MSD by a linear function.
%
% obj.fitMeanMSD computes and fits the weighted mean MSD by a
% straight line y = a * x. The fit is therefore valid only for
% purely diffusive behavior. Fit results are displayed in the
% command window.
%
% obj.fitMeanMSD(clip_factor) does the fit, taking into account
% only the first potion of the average MSD curve specified by
% 'clip_factor' (a double between 0 and 1). If the value
% exceeds 1, then the clip factor is understood to be the
% maximal number of point to take into account in the fit. By
% default, it is set to 0.25.
%
% [fo, gof] = obj.fitMeanMSD(...) returns the fit object and the
% goodness of fit.

if nargin < 2
    clip_factor = 0.25;
end

if ~obj.msd_valid
    obj = obj.computeMSD;
end

ft = fittype('poly1');
mmsd = obj.getMeanMSD;

t = mmsd(:,1);
y = mmsd(:,2);
w = 1./mmsd(:,3);

% Clip data, never take the first one dt = 0
if clip_factor < 1
    t_limit = 2 : round(numel(t) * clip_factor);
else
    t_limit = 1 : min(1+round(clip_factor), numel(t)); %Note: I changed 2 to 1
end
t = t(t_limit);
y = y(t_limit);
w = w(t_limit);

%[fo, gof] = fit(t, y, ft, 'Weights', w); 
[fo, gof] = fit(t, y, ft); %Note: I got rid of weights

str = sprintf([
    'Estimating D through linear fit of the EATA MSD curve (first two delays).\n', ...
    'D = %.3e'], fo.p1/2/obj.n_dim);
disp(str)

if nargout > 0
    varargout{1} = fo;
    if nargout > 1
        varargout{2} = gof;
    end
end

end